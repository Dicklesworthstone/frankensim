//! Typed gesture schedules (music bead
//! `frankensim-music-v8-root-3ez8g.2.3`): gesture is GEOMETRY IN TIME,
//! not MIDI. A performance is a replayable, seedable artifact — "valve 1
//! down at t = 2.10 s over 30 ms; lip pre-stress ramp; sustain pedal at
//! t = 4.0 s; fret the string at 0.32 m at t = 1.2 s" — data the render
//! APIs consume as between-block control deltas, never ad-hoc closures
//! in test files.
//!
//! UNIT DISCIPLINE is type-level: every control target names the ONE
//! [`GestureValue`] variant it accepts (a pressure literally cannot
//! enter a length target — the mismatch refuses at admission, by name).
//! Field names carry their fs-qty dimension in the suffix, the crate's
//! established convention.
//!
//! DETERMINISM: sampling happens on a fixed integer control clock
//! (`control_rate_hz` ticks; no wall time anywhere). Continuous targets
//! ramp linearly over their transition; discrete events (hammer
//! strikes, fret actions, termination swaps) belong to exactly one tick.
//! Canonical bytes + a domain-separated content hash make a schedule a
//! receipt-able artifact.
//!
//! BOUNDARY (binding, from the plan): a valve gesture is an INSERTED
//! LENGTH, a key gesture is a HOLE STATE, a fret gesture is a LENGTH +
//! OBSTACLE — never a mechanism simulation. MIDI import is deliberately
//! absent; note numbers never become omega anywhere in a physics path.

use fs_blake3::{ContentHash, DomainHasher};

/// Schema-versioned hash domain for schedule identities.
pub const GESTURE_SCHEDULE_HASH_DOMAIN: &str = "org.frankensim.fs-scenario.gesture-schedule.v1";
/// Canonical-bytes schema line.
pub const GESTURE_SCHEDULE_SCHEMA: &str = "frankensim-gesture-schedule-v1";

/// What a track controls. Each variant documents the ONE value variant
/// it accepts and the fs-qty dimension of that value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GestureTarget {
    /// Valve crook insertion (value: `LengthM`, dimension L). The length
    /// is the crook-delta receipt's number; 0 = valve up.
    ValveInsertedLength {
        /// Valve index (three bits + tuning slide by convention).
        valve: u8,
    },
    /// Tone-hole opening fraction (value: `Fraction` in [0, 1];
    /// register vents are just more holes).
    ToneHoleSigma {
        /// Hole index from the mouthpiece.
        hole: u8,
    },
    /// Continuous bore-side slide length (value: `LengthM`).
    SlideLength,
    /// String tension trajectory (value: `TensionN`) — retunes and bends
    /// are tension STATE on the string card, never a frequency.
    StringTension {
        /// String index.
        string: u8,
    },
    /// Fret action on a string (value: `Fret` — speaking-length change
    /// AND obstacle engagement; hammer-on/pull-off carry velocity).
    Fret {
        /// String index.
        string: u8,
    },
    /// Bow trajectory (value: `Bow`).
    BowStroke {
        /// String index.
        string: u8,
    },
    /// Piano key strike (value: `StrikeVelocity`; an EVENT).
    HammerStrike {
        /// String/key index.
        string: u8,
    },
    /// Sustain pedal (value: `Fraction`; termination-coupling state,
    /// half-pedal legal).
    SustainPedal,
    /// Una corda (value: `Fraction`; hammer shift / struck-string-count
    /// change on the SAME cards).
    UnaCordaPedal,
    /// Sostenuto hold for one string (value: `Fraction`).
    SostenutoHold {
        /// String index.
        string: u8,
    },
    /// Termination swap by content digest (value: `TerminationSwap`;
    /// mute / hand-in-bell — the bake side is zolja's driver, this is
    /// the control-rate swap event with its D17 fade).
    TerminationSwap,
    /// Blowing / mouth pressure (value: `PressurePa`).
    BlowingPressure,
    /// Lip or fold pre-stress multiplier (value: `Fraction`-like
    /// dimensionless multiplier; > 0).
    PreStress,
    /// Rest aperture of the valve/reed/fold (value: `LengthM`).
    RestAperture,
    /// Flue jet speed (value: `VelocityMPerS`).
    JetSpeed,
    /// Flue jet angle (value: `AngleRad`).
    JetAngle,
}

impl GestureTarget {
    /// Whether `value` is the variant this target accepts (the
    /// type-level unit law).
    #[must_use]
    pub fn accepts(&self, value: &GestureValue) -> bool {
        matches!(
            (self, value),
            (
                GestureTarget::ValveInsertedLength { .. }
                    | GestureTarget::SlideLength
                    | GestureTarget::RestAperture,
                GestureValue::LengthM(_)
            ) | (
                GestureTarget::ToneHoleSigma { .. }
                    | GestureTarget::SustainPedal
                    | GestureTarget::UnaCordaPedal
                    | GestureTarget::SostenutoHold { .. },
                GestureValue::Fraction(_)
            ) | (GestureTarget::PreStress, GestureValue::Multiplier(_))
                | (
                    GestureTarget::StringTension { .. },
                    GestureValue::TensionN(_)
                )
                | (GestureTarget::Fret { .. }, GestureValue::Fret { .. })
                | (GestureTarget::BowStroke { .. }, GestureValue::Bow { .. })
                | (
                    GestureTarget::HammerStrike { .. },
                    GestureValue::StrikeVelocity { .. }
                )
                | (
                    GestureTarget::TerminationSwap,
                    GestureValue::TerminationSwap { .. }
                )
                | (GestureTarget::BlowingPressure, GestureValue::PressurePa(_))
                | (GestureTarget::JetSpeed, GestureValue::VelocityMPerS(_))
                | (GestureTarget::JetAngle, GestureValue::AngleRad(_))
        )
    }

    /// Whether this target is EVENT-like (belongs to one tick) rather
    /// than a ramped continuous quantity.
    #[must_use]
    pub fn is_event(&self) -> bool {
        matches!(
            self,
            GestureTarget::HammerStrike { .. }
                | GestureTarget::Fret { .. }
                | GestureTarget::TerminationSwap
        )
    }
}

/// A typed control value. Dimensions are in the variant, not a tag.
#[derive(Debug, Clone, PartialEq)]
pub enum GestureValue {
    /// Length [m].
    LengthM(f64),
    /// Fraction in [0, 1].
    Fraction(f64),
    /// Dimensionless positive multiplier.
    Multiplier(f64),
    /// Tension [N].
    TensionN(f64),
    /// Pressure [Pa].
    PressurePa(f64),
    /// Velocity [m/s].
    VelocityMPerS(f64),
    /// Angle [rad].
    AngleRad(f64),
    /// Fret action: length change AND obstacle engagement.
    Fret {
        /// Engaged (hammer-on / press) or released (pull-off).
        engaged: bool,
        /// Fret position from the nut [m] (the new speaking length).
        position_m: f64,
        /// Obstacle height under the string [m].
        height_m: f64,
        /// Action velocity [m/s] (hammer-on/pull-off strength).
        velocity_m_per_s: f64,
    },
    /// Bow state.
    Bow {
        /// Bow velocity [m/s].
        velocity_m_per_s: f64,
        /// Normal force [N].
        normal_force_n: f64,
        /// Bowing station as a fraction of speaking length in (0, 1).
        station: f64,
    },
    /// Hammer strike velocity [m/s].
    StrikeVelocity {
        /// Impact velocity [m/s].
        velocity_m_per_s: f64,
    },
    /// Termination table swap by content digest with a D17 fade.
    TerminationSwap {
        /// Hex digest of the target tabulated load.
        digest_hex: String,
        /// Crossfade duration [s].
        fade_s: f64,
    },
}

impl GestureValue {
    fn finite(&self) -> bool {
        match self {
            GestureValue::LengthM(v)
            | GestureValue::Fraction(v)
            | GestureValue::Multiplier(v)
            | GestureValue::TensionN(v)
            | GestureValue::PressurePa(v)
            | GestureValue::VelocityMPerS(v)
            | GestureValue::AngleRad(v) => v.is_finite(),
            GestureValue::Fret {
                position_m,
                height_m,
                velocity_m_per_s,
                ..
            } => position_m.is_finite() && height_m.is_finite() && velocity_m_per_s.is_finite(),
            GestureValue::Bow {
                velocity_m_per_s,
                normal_force_n,
                station,
            } => velocity_m_per_s.is_finite() && normal_force_n.is_finite() && station.is_finite(),
            GestureValue::StrikeVelocity { velocity_m_per_s } => velocity_m_per_s.is_finite(),
            GestureValue::TerminationSwap { fade_s, .. } => fade_s.is_finite(),
        }
    }

    fn in_range(&self) -> bool {
        match self {
            GestureValue::Fraction(v) => (0.0..=1.0).contains(v),
            GestureValue::Multiplier(v) => *v > 0.0,
            GestureValue::LengthM(v) | GestureValue::TensionN(v) | GestureValue::PressurePa(v) => {
                *v >= 0.0
            }
            GestureValue::Bow { station, .. } => *station > 0.0 && *station < 1.0,
            GestureValue::TerminationSwap { fade_s, digest_hex } => {
                *fade_s >= 0.0 && !digest_hex.is_empty()
            }
            _ => true,
        }
    }

    /// Scalar payload for continuous interpolation (events return None).
    fn scalar(&self) -> Option<f64> {
        match self {
            GestureValue::LengthM(v)
            | GestureValue::Fraction(v)
            | GestureValue::Multiplier(v)
            | GestureValue::TensionN(v)
            | GestureValue::PressurePa(v)
            | GestureValue::VelocityMPerS(v)
            | GestureValue::AngleRad(v) => Some(*v),
            _ => None,
        }
    }
}

/// One scheduled change on a track.
#[derive(Debug, Clone, PartialEq)]
pub struct GestureEvent {
    /// Start time [s] on the control clock (>= 0, strictly increasing
    /// per track).
    pub time_s: f64,
    /// Ramp duration [s] to reach the value (0 = step; must be 0 for
    /// event-like targets).
    pub transition_s: f64,
    /// The typed value.
    pub value: GestureValue,
}

/// One control track: a named target plus its ordered events.
#[derive(Debug, Clone, PartialEq)]
pub struct GestureTrack {
    /// Unique id the consumer binds to (e.g. `valve-1`, `blow`).
    pub id: String,
    /// The typed target.
    pub target: GestureTarget,
    /// Initial value (before the first event).
    pub initial: GestureValue,
    /// Ordered events.
    pub events: Vec<GestureEvent>,
}

/// A complete replayable performance schedule.
#[derive(Debug, Clone, PartialEq)]
pub struct GestureSchedule {
    /// Control clock [Hz] (fixed; no wall time anywhere).
    pub control_rate_hz: u32,
    tracks: Vec<GestureTrack>,
}

/// Typed refusals — nothing silently no-ops.
#[derive(Debug, Clone, PartialEq)]
pub enum GestureError {
    /// A structural parameter is unusable.
    Invalid {
        /// Diagnosis.
        what: &'static str,
    },
    /// A value variant does not match its target (the unit law).
    UnitMismatch {
        /// The offending track id.
        track: String,
    },
    /// A value is outside its legal range (sigma, station, ...).
    OutOfRange {
        /// The offending track id.
        track: String,
    },
    /// Track times are not strictly increasing.
    NonMonotoneTime {
        /// The offending track id.
        track: String,
    },
    /// A consumer asked for a control id the schedule does not carry —
    /// refused BY NAME, never a silent no-op (the falsifier's arm).
    UnknownControlId {
        /// The requested id.
        requested: String,
    },
    /// Canonical bytes failed decode at a named line.
    Decode {
        /// 1-based line.
        line: usize,
    },
}

impl core::fmt::Display for GestureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            GestureError::Invalid { what } => write!(f, "FS-SCENARIO-GESTURE: {what}"),
            GestureError::UnitMismatch { track } => write!(
                f,
                "FS-SCENARIO-GESTURE-UNIT: track {track:?} value variant does not match its \
                 target (a pressure cannot enter a length control)"
            ),
            GestureError::OutOfRange { track } => {
                write!(
                    f,
                    "FS-SCENARIO-GESTURE-RANGE: track {track:?} value out of range"
                )
            }
            GestureError::NonMonotoneTime { track } => write!(
                f,
                "FS-SCENARIO-GESTURE-TIME: track {track:?} event times must strictly increase"
            ),
            GestureError::UnknownControlId { requested } => write!(
                f,
                "FS-SCENARIO-GESTURE-ID: no track named {requested:?} (refused by name, never \
                 a silent no-op)"
            ),
            GestureError::Decode { line } => {
                write!(f, "FS-SCENARIO-GESTURE-DECODE: line {line}")
            }
        }
    }
}

impl core::error::Error for GestureError {}

impl GestureSchedule {
    /// Admit a schedule: per-track unit law, ranges, monotone times,
    /// unique ids, event-like targets step-only.
    ///
    /// # Errors
    /// [`GestureError`] naming the violated law and track.
    pub fn try_new(
        control_rate_hz: u32,
        tracks: Vec<GestureTrack>,
    ) -> Result<GestureSchedule, GestureError> {
        if control_rate_hz == 0 {
            return Err(GestureError::Invalid {
                what: "control rate must be positive",
            });
        }
        let mut seen = std::collections::BTreeSet::new();
        for track in &tracks {
            if track.id.trim().is_empty() {
                return Err(GestureError::Invalid {
                    what: "track ids must be non-empty",
                });
            }
            if !seen.insert(track.id.clone()) {
                return Err(GestureError::Invalid {
                    what: "track ids must be unique",
                });
            }
            let check_value = |value: &GestureValue| -> Result<(), GestureError> {
                if !track.target.accepts(value) {
                    return Err(GestureError::UnitMismatch {
                        track: track.id.clone(),
                    });
                }
                if !value.finite() || !value.in_range() {
                    return Err(GestureError::OutOfRange {
                        track: track.id.clone(),
                    });
                }
                Ok(())
            };
            check_value(&track.initial)?;
            let mut last = -1.0f64;
            for event in &track.events {
                if !(event.time_s.is_finite() && event.time_s >= 0.0 && event.time_s > last) {
                    return Err(GestureError::NonMonotoneTime {
                        track: track.id.clone(),
                    });
                }
                last = event.time_s;
                if !(event.transition_s.is_finite() && event.transition_s >= 0.0) {
                    return Err(GestureError::Invalid {
                        what: "transitions must be finite non-negative",
                    });
                }
                if track.target.is_event() && event.transition_s != 0.0 {
                    return Err(GestureError::Invalid {
                        what: "event-like targets take step events only",
                    });
                }
                check_value(&event.value)?;
            }
        }
        Ok(GestureSchedule {
            control_rate_hz,
            tracks,
        })
    }

    /// The track list (read-only; admission owns mutation).
    #[must_use]
    pub fn tracks(&self) -> &[GestureTrack] {
        &self.tracks
    }

    /// Sample a CONTINUOUS track at control tick `tick` (linear ramp
    /// over each event's transition; deterministic pure function of the
    /// integer tick).
    ///
    /// # Errors
    /// [`GestureError::UnknownControlId`]; `Invalid` for event tracks.
    pub fn sample(&self, id: &str, tick: u64) -> Result<f64, GestureError> {
        let track = self.tracks.iter().find(|t| t.id == id).ok_or_else(|| {
            GestureError::UnknownControlId {
                requested: id.to_string(),
            }
        })?;
        if track.target.is_event() {
            return Err(GestureError::Invalid {
                what: "event tracks are read with events_at, not sample",
            });
        }
        let t = tick as f64 / f64::from(self.control_rate_hz);
        let mut value = track.initial.scalar().expect("continuous by admission");
        for event in &track.events {
            let target = event.value.scalar().expect("continuous by admission");
            if t < event.time_s {
                break;
            }
            if event.transition_s > 0.0 && t < event.time_s + event.transition_s {
                let f = (t - event.time_s) / event.transition_s;
                value += (target - value) * f;
                break;
            }
            value = target;
        }
        Ok(value)
    }

    /// The EVENTS of an event-like track that belong to control tick
    /// `tick` (each event belongs to exactly one tick: `floor(time *
    /// rate) == tick`).
    ///
    /// # Errors
    /// [`GestureError::UnknownControlId`]; `Invalid` for continuous
    /// tracks.
    pub fn events_at(&self, id: &str, tick: u64) -> Result<Vec<&GestureValue>, GestureError> {
        let track = self.tracks.iter().find(|t| t.id == id).ok_or_else(|| {
            GestureError::UnknownControlId {
                requested: id.to_string(),
            }
        })?;
        if !track.target.is_event() {
            return Err(GestureError::Invalid {
                what: "continuous tracks are read with sample, not events_at",
            });
        }
        let rate = f64::from(self.control_rate_hz);
        Ok(track
            .events
            .iter()
            .filter(|e| (e.time_s * rate).floor() as u64 == tick)
            .map(|e| &e.value)
            .collect())
    }

    /// Canonical line-oriented bytes (schema first, fixed field order,
    /// `{:e}` floats; no wall time, no commit stamps).
    #[must_use]
    pub fn to_canonical_bytes(&self) -> Vec<u8> {
        let mut s = String::new();
        s.push_str(GESTURE_SCHEDULE_SCHEMA);
        s.push('\n');
        s.push_str(&format!("control_rate_hz\t{}\n", self.control_rate_hz));
        s.push_str(&format!("tracks\t{}\n", self.tracks.len()));
        for track in &self.tracks {
            s.push_str(&format!(
                "track\t{}\t{}\n",
                track.id,
                encode_target(&track.target)
            ));
            s.push_str(&format!("initial\t{}\n", encode_value(&track.initial)));
            s.push_str(&format!("events\t{}\n", track.events.len()));
            for e in &track.events {
                s.push_str(&format!(
                    "event\t{:e}\t{:e}\t{}\n",
                    e.time_s,
                    e.transition_s,
                    encode_value(&e.value)
                ));
            }
        }
        s.into_bytes()
    }

    /// Domain-separated content hash of the canonical bytes.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        let mut h = DomainHasher::new(GESTURE_SCHEDULE_HASH_DOMAIN);
        h.update(&self.to_canonical_bytes());
        h.finalize()
    }

    /// Decode + re-admit canonical bytes (all admission laws re-run).
    ///
    /// # Errors
    /// [`GestureError::Decode`] at the first bad line; admission errors
    /// as in [`GestureSchedule::try_new`].
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<GestureSchedule, GestureError> {
        let text = core::str::from_utf8(bytes).map_err(|_| GestureError::Decode { line: 1 })?;
        let lines: Vec<&str> = text.lines().collect();
        let mut at = 0usize;
        let mut next = |at: &mut usize| -> Result<&str, GestureError> {
            let line = lines
                .get(*at)
                .copied()
                .ok_or(GestureError::Decode { line: *at + 1 })?;
            *at += 1;
            Ok(line)
        };
        if next(&mut at)? != GESTURE_SCHEDULE_SCHEMA {
            return Err(GestureError::Decode { line: 1 });
        }
        let rate_line = next(&mut at)?;
        let control_rate_hz = rate_line
            .strip_prefix("control_rate_hz\t")
            .and_then(|v| v.parse::<u32>().ok())
            .ok_or(GestureError::Decode { line: at })?;
        let n_tracks = next(&mut at)?
            .strip_prefix("tracks\t")
            .and_then(|v| v.parse::<usize>().ok())
            .ok_or(GestureError::Decode { line: at })?;
        let mut tracks = Vec::with_capacity(n_tracks);
        for _ in 0..n_tracks {
            let header = next(&mut at)?;
            let mut cols = header.split('\t');
            if cols.next() != Some("track") {
                return Err(GestureError::Decode { line: at });
            }
            let id = cols
                .next()
                .ok_or(GestureError::Decode { line: at })?
                .to_string();
            let target = decode_target(&cols.collect::<Vec<_>>().join("\t"))
                .ok_or(GestureError::Decode { line: at })?;
            let initial_line = next(&mut at)?;
            let initial = initial_line
                .strip_prefix("initial\t")
                .and_then(decode_value)
                .ok_or(GestureError::Decode { line: at })?;
            let n_events = next(&mut at)?
                .strip_prefix("events\t")
                .and_then(|v| v.parse::<usize>().ok())
                .ok_or(GestureError::Decode { line: at })?;
            let mut events = Vec::with_capacity(n_events);
            for _ in 0..n_events {
                let line = next(&mut at)?;
                let rest = line
                    .strip_prefix("event\t")
                    .ok_or(GestureError::Decode { line: at })?;
                let mut cols = rest.splitn(3, '\t');
                let time_s = cols
                    .next()
                    .and_then(|v| v.parse::<f64>().ok())
                    .ok_or(GestureError::Decode { line: at })?;
                let transition_s = cols
                    .next()
                    .and_then(|v| v.parse::<f64>().ok())
                    .ok_or(GestureError::Decode { line: at })?;
                let value = cols
                    .next()
                    .and_then(decode_value)
                    .ok_or(GestureError::Decode { line: at })?;
                events.push(GestureEvent {
                    time_s,
                    transition_s,
                    value,
                });
            }
            tracks.push(GestureTrack {
                id,
                target,
                initial,
                events,
            });
        }
        GestureSchedule::try_new(control_rate_hz, tracks)
    }
}

fn encode_target(t: &GestureTarget) -> String {
    match t {
        GestureTarget::ValveInsertedLength { valve } => format!("valve-length\t{valve}"),
        GestureTarget::ToneHoleSigma { hole } => format!("hole-sigma\t{hole}"),
        GestureTarget::SlideLength => "slide-length".to_string(),
        GestureTarget::StringTension { string } => format!("string-tension\t{string}"),
        GestureTarget::Fret { string } => format!("fret\t{string}"),
        GestureTarget::BowStroke { string } => format!("bow\t{string}"),
        GestureTarget::HammerStrike { string } => format!("hammer\t{string}"),
        GestureTarget::SustainPedal => "sustain".to_string(),
        GestureTarget::UnaCordaPedal => "una-corda".to_string(),
        GestureTarget::SostenutoHold { string } => format!("sostenuto\t{string}"),
        GestureTarget::TerminationSwap => "termination-swap".to_string(),
        GestureTarget::BlowingPressure => "blowing-pressure".to_string(),
        GestureTarget::PreStress => "pre-stress".to_string(),
        GestureTarget::RestAperture => "rest-aperture".to_string(),
        GestureTarget::JetSpeed => "jet-speed".to_string(),
        GestureTarget::JetAngle => "jet-angle".to_string(),
    }
}

fn decode_target(s: &str) -> Option<GestureTarget> {
    let mut cols = s.split('\t');
    let kind = cols.next()?;
    let idx = |cols: &mut core::str::Split<'_, char>| cols.next()?.parse::<u8>().ok();
    Some(match kind {
        "valve-length" => GestureTarget::ValveInsertedLength {
            valve: idx(&mut cols)?,
        },
        "hole-sigma" => GestureTarget::ToneHoleSigma {
            hole: idx(&mut cols)?,
        },
        "slide-length" => GestureTarget::SlideLength,
        "string-tension" => GestureTarget::StringTension {
            string: idx(&mut cols)?,
        },
        "fret" => GestureTarget::Fret {
            string: idx(&mut cols)?,
        },
        "bow" => GestureTarget::BowStroke {
            string: idx(&mut cols)?,
        },
        "hammer" => GestureTarget::HammerStrike {
            string: idx(&mut cols)?,
        },
        "sustain" => GestureTarget::SustainPedal,
        "una-corda" => GestureTarget::UnaCordaPedal,
        "sostenuto" => GestureTarget::SostenutoHold {
            string: idx(&mut cols)?,
        },
        "termination-swap" => GestureTarget::TerminationSwap,
        "blowing-pressure" => GestureTarget::BlowingPressure,
        "pre-stress" => GestureTarget::PreStress,
        "rest-aperture" => GestureTarget::RestAperture,
        "jet-speed" => GestureTarget::JetSpeed,
        "jet-angle" => GestureTarget::JetAngle,
        _ => return None,
    })
}

fn encode_value(v: &GestureValue) -> String {
    match v {
        GestureValue::LengthM(x) => format!("length-m\t{x:e}"),
        GestureValue::Fraction(x) => format!("fraction\t{x:e}"),
        GestureValue::Multiplier(x) => format!("multiplier\t{x:e}"),
        GestureValue::TensionN(x) => format!("tension-n\t{x:e}"),
        GestureValue::PressurePa(x) => format!("pressure-pa\t{x:e}"),
        GestureValue::VelocityMPerS(x) => format!("velocity-m-s\t{x:e}"),
        GestureValue::AngleRad(x) => format!("angle-rad\t{x:e}"),
        GestureValue::Fret {
            engaged,
            position_m,
            height_m,
            velocity_m_per_s,
        } => format!(
            "fret\t{}\t{position_m:e}\t{height_m:e}\t{velocity_m_per_s:e}",
            u8::from(*engaged)
        ),
        GestureValue::Bow {
            velocity_m_per_s,
            normal_force_n,
            station,
        } => format!("bow\t{velocity_m_per_s:e}\t{normal_force_n:e}\t{station:e}"),
        GestureValue::StrikeVelocity { velocity_m_per_s } => {
            format!("strike\t{velocity_m_per_s:e}")
        }
        GestureValue::TerminationSwap { digest_hex, fade_s } => {
            format!("termination-swap\t{digest_hex}\t{fade_s:e}")
        }
    }
}

fn decode_value(s: &str) -> Option<GestureValue> {
    let mut cols = s.split('\t');
    let kind = cols.next()?;
    let num = |cols: &mut core::str::Split<'_, char>| cols.next()?.parse::<f64>().ok();
    Some(match kind {
        "length-m" => GestureValue::LengthM(num(&mut cols)?),
        "fraction" => GestureValue::Fraction(num(&mut cols)?),
        "multiplier" => GestureValue::Multiplier(num(&mut cols)?),
        "tension-n" => GestureValue::TensionN(num(&mut cols)?),
        "pressure-pa" => GestureValue::PressurePa(num(&mut cols)?),
        "velocity-m-s" => GestureValue::VelocityMPerS(num(&mut cols)?),
        "angle-rad" => GestureValue::AngleRad(num(&mut cols)?),
        "fret" => GestureValue::Fret {
            engaged: cols.next()? == "1",
            position_m: num(&mut cols)?,
            height_m: num(&mut cols)?,
            velocity_m_per_s: num(&mut cols)?,
        },
        "bow" => GestureValue::Bow {
            velocity_m_per_s: num(&mut cols)?,
            normal_force_n: num(&mut cols)?,
            station: num(&mut cols)?,
        },
        "strike" => GestureValue::StrikeVelocity {
            velocity_m_per_s: num(&mut cols)?,
        },
        "termination-swap" => GestureValue::TerminationSwap {
            digest_hex: cols.next()?.to_string(),
            fade_s: num(&mut cols)?,
        },
        _ => return None,
    })
}

#[cfg(test)]
mod gesture_tests {
    use super::*;

    fn verdict(case: &str, pass: bool, detail: &str) {
        println!(
            "{{\"suite\":\"fs-scenario\",\"case\":\"{case}\",\"verdict\":\"{}\",\"detail\":\"{detail}\"}}",
            if pass { "pass" } else { "fail" }
        );
        assert!(pass, "case {case}: {detail}");
    }

    fn performance() -> GestureSchedule {
        GestureSchedule::try_new(
            200,
            vec![
                GestureTrack {
                    id: "blow".to_string(),
                    target: GestureTarget::BlowingPressure,
                    initial: GestureValue::PressurePa(0.0),
                    events: vec![
                        GestureEvent {
                            time_s: 0.1,
                            transition_s: 0.05,
                            value: GestureValue::PressurePa(3000.0),
                        },
                        GestureEvent {
                            time_s: 1.5,
                            transition_s: 0.2,
                            value: GestureValue::PressurePa(1000.0),
                        },
                    ],
                },
                GestureTrack {
                    id: "valve-1".to_string(),
                    target: GestureTarget::ValveInsertedLength { valve: 1 },
                    initial: GestureValue::LengthM(0.0),
                    events: vec![GestureEvent {
                        time_s: 2.10,
                        transition_s: 0.03,
                        value: GestureValue::LengthM(0.16),
                    }],
                },
                GestureTrack {
                    id: "hole-3".to_string(),
                    target: GestureTarget::ToneHoleSigma { hole: 3 },
                    initial: GestureValue::Fraction(1.0),
                    events: vec![GestureEvent {
                        time_s: 0.8,
                        transition_s: 0.01,
                        value: GestureValue::Fraction(0.0),
                    }],
                },
                GestureTrack {
                    id: "fret-2".to_string(),
                    target: GestureTarget::Fret { string: 2 },
                    initial: GestureValue::Fret {
                        engaged: false,
                        position_m: 0.0,
                        height_m: 0.0,
                        velocity_m_per_s: 0.0,
                    },
                    events: vec![GestureEvent {
                        time_s: 1.2,
                        transition_s: 0.0,
                        value: GestureValue::Fret {
                            engaged: true,
                            position_m: 0.32,
                            height_m: 1.5e-3,
                            velocity_m_per_s: 0.8,
                        },
                    }],
                },
                GestureTrack {
                    id: "strike-40".to_string(),
                    target: GestureTarget::HammerStrike { string: 40 },
                    initial: GestureValue::StrikeVelocity {
                        velocity_m_per_s: 0.0,
                    },
                    events: vec![GestureEvent {
                        time_s: 3.0,
                        transition_s: 0.0,
                        value: GestureValue::StrikeVelocity {
                            velocity_m_per_s: 2.5,
                        },
                    }],
                },
                GestureTrack {
                    id: "mute".to_string(),
                    target: GestureTarget::TerminationSwap,
                    initial: GestureValue::TerminationSwap {
                        digest_hex: "open".to_string(),
                        fade_s: 0.0,
                    },
                    events: vec![GestureEvent {
                        time_s: 4.0,
                        transition_s: 0.0,
                        value: GestureValue::TerminationSwap {
                            digest_hex: "ab12cd34".to_string(),
                            fade_s: 0.004,
                        },
                    }],
                },
            ],
        )
        .expect("performance admits")
    }

    #[test]
    fn gs_001_deterministic_sampling_and_ramps() {
        let s = performance();
        // Before, mid-ramp, and settled samples are pure functions of the
        // integer tick.
        let before = s.sample("blow", 10).expect("t=0.05");
        let mid = s.sample("blow", 25).expect("t=0.125");
        let after = s.sample("blow", 60).expect("t=0.30");
        let down_mid = s.sample("blow", 320).expect("t=1.6");
        let bitwise = (0..1000u64)
            .all(|k| s.sample("blow", k).expect("a") == s.sample("blow", k).expect("b"));
        // Hole closes 1 -> 0 at 0.8 s.
        let open = s.sample("hole-3", 100).expect("t=0.5");
        let closed = s.sample("hole-3", 200).expect("t=1.0");
        // Events belong to exactly one tick.
        let strike_tick = (3.0f64 * 200.0).floor() as u64;
        let hits = s.events_at("strike-40", strike_tick).expect("events");
        let misses = s.events_at("strike-40", strike_tick + 1).expect("events");
        let fret_hits = s.events_at("fret-2", 240).expect("fret tick");
        let pass = before == 0.0
            && mid > 1000.0
            && mid < 2000.0
            && after == 3000.0
            && down_mid > 1000.0
            && down_mid < 3000.0
            && bitwise
            && open == 1.0
            && closed == 0.0
            && hits.len() == 1
            && misses.is_empty()
            && fret_hits.len() == 1;
        verdict(
            "gs-001-deterministic-sampling",
            pass,
            &format!(
                "blow 0/{mid:.0}/3000, down-ramp {down_mid:.0}, bitwise {bitwise}, hole \
                 {open}->{closed}, strike tick {strike_tick} hits {}, fret hits {}",
                hits.len(),
                fret_hits.len()
            ),
        );
    }

    #[test]
    fn gs_002_round_trip_and_content_hash() {
        let s = performance();
        let bytes = s.to_canonical_bytes();
        let back = GestureSchedule::from_canonical_bytes(&bytes).expect("round trip");
        let identical = back == s;
        let hash_stable = back.content_hash() == s.content_hash();
        // A one-byte edit changes the hash and/or refuses.
        let mut tampered = bytes.clone();
        let pos = bytes
            .windows(4)
            .position(|w| w == b"3e3\n" || w == b"3e0\n")
            .unwrap_or(40);
        tampered[pos] ^= 1;
        let tampered_differs = match GestureSchedule::from_canonical_bytes(&tampered) {
            Ok(t) => t.content_hash() != s.content_hash(),
            Err(_) => true,
        };
        let pass = identical && hash_stable && tampered_differs;
        verdict(
            "gs-002-round-trip",
            pass,
            &format!(
                "round-trip identical {identical}, hash stable {hash_stable}, tamper \
                 detected {tampered_differs}, hash {}",
                s.content_hash().to_hex()
            ),
        );
    }

    #[test]
    fn gs_003_refusals_by_name() {
        // Unit mismatch: a pressure into a length target.
        let unit = GestureSchedule::try_new(
            100,
            vec![GestureTrack {
                id: "v".to_string(),
                target: GestureTarget::ValveInsertedLength { valve: 1 },
                initial: GestureValue::PressurePa(100.0),
                events: vec![],
            }],
        );
        // Out-of-range sigma.
        let range = GestureSchedule::try_new(
            100,
            vec![GestureTrack {
                id: "h".to_string(),
                target: GestureTarget::ToneHoleSigma { hole: 0 },
                initial: GestureValue::Fraction(1.4),
                events: vec![],
            }],
        );
        // Non-monotone time.
        let time = GestureSchedule::try_new(
            100,
            vec![GestureTrack {
                id: "b".to_string(),
                target: GestureTarget::BlowingPressure,
                initial: GestureValue::PressurePa(0.0),
                events: vec![
                    GestureEvent {
                        time_s: 1.0,
                        transition_s: 0.0,
                        value: GestureValue::PressurePa(1.0),
                    },
                    GestureEvent {
                        time_s: 0.5,
                        transition_s: 0.0,
                        value: GestureValue::PressurePa(2.0),
                    },
                ],
            }],
        );
        // Ramped event-target refuses.
        let ramped_event = GestureSchedule::try_new(
            100,
            vec![GestureTrack {
                id: "s".to_string(),
                target: GestureTarget::HammerStrike { string: 0 },
                initial: GestureValue::StrikeVelocity {
                    velocity_m_per_s: 0.0,
                },
                events: vec![GestureEvent {
                    time_s: 1.0,
                    transition_s: 0.5,
                    value: GestureValue::StrikeVelocity {
                        velocity_m_per_s: 1.0,
                    },
                }],
            }],
        );
        // THE FALSIFIER: an unknown control id refuses BY NAME, never a
        // silent no-op.
        let s = performance();
        let unknown = s.sample("does-not-exist", 0);
        let pass = matches!(unit, Err(GestureError::UnitMismatch { .. }))
            && matches!(range, Err(GestureError::OutOfRange { .. }))
            && matches!(time, Err(GestureError::NonMonotoneTime { .. }))
            && matches!(ramped_event, Err(GestureError::Invalid { .. }))
            && matches!(unknown, Err(GestureError::UnknownControlId { .. }));
        verdict(
            "gs-003-refusals",
            pass,
            &format!(
                "unit {} range {} time {} ramped-event {} unknown-id {}",
                unit.is_err(),
                range.is_err(),
                time.is_err(),
                ramped_event.is_err(),
                unknown.is_err()
            ),
        );
    }

    #[test]
    fn gs_004_fret_is_length_plus_obstacle() {
        // The fret event carries BOTH sides of the plan's row: the new
        // speaking length AND the obstacle engagement — applied to a
        // string description, both are logged.
        let s = performance();
        let fret = s.events_at("fret-2", 240).expect("fret tick");
        let GestureValue::Fret {
            engaged,
            position_m,
            height_m,
            velocity_m_per_s,
        } = fret[0]
        else {
            panic!("fret event variant");
        };
        let string = crate::PrestressedString {
            length_m: 0.65,
            tension_n: 60.0,
            lin_density_kg_m: 6.0e-4,
            axial_stiffness_n: 0.0,
            width_m: 1.0e-3,
            n_modes: 8,
            damping_ratio: 1.0e-3,
            rayleigh: None,
            bending_stiffness_n_m2: 0.0,
            polarization_detune: 0.0,
            moving_end: false,
        };
        // Length side: the fret shortens the SPEAKING length.
        let speaking = if *engaged {
            *position_m
        } else {
            string.length_m
        };
        // Obstacle side: engagement parameters for the fs-dcontact lay.
        let pass =
            *engaged && speaking < string.length_m && *height_m > 0.0 && *velocity_m_per_s > 0.0;
        verdict(
            "gs-004-fret-both-sides",
            pass,
            &format!(
                "fret engaged {engaged}: speaking length {:.3} -> {speaking:.3} m AND \
                 obstacle (height {height_m:.1e} m, action velocity {velocity_m_per_s:.1} \
                 m/s) — length + obstacle, never a mechanism simulation",
                string.length_m
            ),
        );
    }
}

//! File-driven compliant stops and normal trajectories in the original material/bow/contact performer.
//! This profile describes authored contact laws, never inferred tissue or fret
//! properties. Existing source-derived centreline and receiver shapes own force
//! projection; the common contact, preload and time solvers own all mechanics.

use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::{
    MultiBowStringBank, stops::{MaterialStringStop, StringStopSide},
};
use fs_couple::render::schedule::force::coupled::contact::ModalContactConfig;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::obstacles::FrictionObstacleConfig;
use fs_dcontact::Obstacle;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::contact_load::BowContactLoad;
use fs_couple::bowed_string::runtime::reciprocal::bank::polarization::TransverseDirection;

pub const MAX_BYTES: usize = 16_384;
const MAX_STOPS: usize = 32;
const MAX_MOTION_KNOTS: usize = 128;
const MAX_MOTION_COORDINATES: usize = MAX_STOPS * MAX_MOTION_KNOTS;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::motion::{
    StopMotionProgram, StopMotionKnot, ObstacleMotionLimits,
};
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::stations::{
    StopStationProgram, StopStationKnot, StopStationMotionLimits,
};

#[derive(Clone, Debug)]
struct Entry {
    name: String,
    location: StringBowLocation,
    side: StringStopSide,
    // None is an explicitly immovable support. Some selects the original body
    // footprint of a named string's sole bridge, NOT that string's displacement.
    receiver_bridge: Option<usize>,
    law: Obstacle,
    limits: ModalContactConfig,
}
#[derive(Clone, Debug)]
pub struct Profile {
    pub bytes: Vec<u8>,
    entries: Vec<Entry>,
    motion: Option<StopMotionProgram>,
    station_motion: Option<StopStationProgram>,
    bow_contacts: Vec<BowContactLoad>,
}

fn scalar(word: &str) -> Result<f64> {
    let value: f64 = word.parse()?;
    if !value.is_finite() { return Err("stop scalar must be finite".into()); }
    Ok(value)
}
impl Profile {
    pub fn parse(bytes: &[u8]) -> Result<Self> {
        if bytes.len() > MAX_BYTES { return Err("stop profile exceeds 16 KiB".into()); }
        let text = std::str::from_utf8(bytes)?;
        let mut lines = text.lines();
        if lines.next() != Some("frankensim-string-stops-v1") {
            return Err("unsupported string-stop profile version".into());
        }
        let mut entries: Vec<Entry> = Vec::new();
        let mut motion = None;
        let mut station_motion = None;
        let mut bow_contacts: Vec<BowContactLoad> = Vec::new();
        while let Some(line) = lines.next() {
            let fields: Vec<_> = line.split_ascii_whitespace().collect();
            if fields.first() == Some(&"motion") {
                if motion.is_some() { return Err("normal motion can be declared only once".into()); }
                motion = Some(parse_motion(&fields, &mut lines, entries.len())?);
                continue;
            }
            if fields.first() == Some(&"station-motion") {
                if station_motion.is_some() { return Err("station motion can be declared only once".into()); }
                station_motion = Some(parse_station_motion(&fields, &mut lines, entries.len())?);
                continue;
            }
            if motion.is_some() || station_motion.is_some() {
                return Err("stop and bow-contact declarations must precede the geometry motion blocks".into());
            }
            if fields.first() == Some(&"bow-contact") {
                if fields.len() != 3 { return Err("expected bow-contact BOW_INDEX STOP_ID".into()); }
                let bow: usize = fields[1].parse()?;
                let stop = entries.iter().position(|e| e.name == fields[2])
                    .ok_or("bow-contact must name a preceding stop declaration")?;
                if bow >= 16 || bow_contacts.iter().any(|b| b.bow == bow || b.stop == stop) {
                    return Err("bow contacts require distinct actual bows and distinct normal contacts".into());
                }
                bow_contacts.push(BowContactLoad { bow, stop });
                continue;
            }
            if entries.len() == MAX_STOPS { return Err("stop profile exceeds 32 contacts".into()); }
            if fields.len() < 5 || fields[0] != "stop" {
                return Err("expected stop ID STRING STATION SIDE RECEIVER GAP K ALPHA LOSS MAX_FORCE MAX_PENETRATION".into());
            }
            let name = fields[1];
            if name.len() > 64 || !name.bytes().all(|b| b.is_ascii_alphanumeric() || b"_-".contains(&b))
                || entries.iter().any(|e| e.name == name) {
                return Err("stop identities must be unique short ASCII tokens".into());
            }
            let location = StringBowLocation { string: fields[2].parse()?, station_fraction: scalar(fields[3])? };
            let (side, values) = match fields[4] {
                "positive" if fields.len() == 12 => (StringStopSide::Positive, &fields[5..]),
                "negative" if fields.len() == 12 => (StringStopSide::Negative, &fields[5..]),
                "normal" if fields.len() == 14 => (StringStopSide::Along(
                    TransverseDirection::try_new([scalar(fields[5])?, scalar(fields[6])?])?), &fields[7..]),
                _ => return Err("stop side must be positive, negative or normal NX NY with a complete record".into()),
            };
            let receiver_bridge = if values[0] == "fixed" { None } else {
                Some(values[0].strip_prefix("bridge:").ok_or("receiver must be fixed or bridge:STRING")?.parse::<usize>()?)
            };
            let gap = scalar(values[1])?;
            let stiffness = scalar(values[2])?;
            let alpha = scalar(values[3])?;
            let loss = scalar(values[4])?;
            let max_force = scalar(values[5])?;
            let max_penetration = scalar(values[6])?;
            if !(0.0 < location.station_fraction && location.station_fraction < 1.0)
                || stiffness <= 0.0 || alpha < 1.0 || loss < 0.0 || max_force <= 0.0 || max_penetration <= 0.0 {
                return Err("stop geometry, law and physical limits are invalid".into());
            }
            // Keep the original point law's units: k is N/m^alpha; unit weight
            // adds no fictitious patch area. Gap is unloaded surface clearance.
            let law = Obstacle::new(vec![-1.0], 1, 1, vec![gap], vec![1.0], stiffness, alpha,
                format!("authored string-stop profile v1: {name}; not measured contact data"))?
                .with_internal_loss(loss)?;
            entries.push(Entry { name: name.into(), location, side, receiver_bridge, law,
                limits: ModalContactConfig { max_iterations: 96, maximum_force_n: max_force,
                    maximum_penetration_m: max_penetration, force_absolute_tolerance_n: 1e-11,
                    force_relative_tolerance: 1e-9 } });
        }
        if entries.is_empty() { return Err("stop profile must declare at least one contact".into()); }
        Ok(Self { bytes: bytes.to_vec(), entries, motion, station_motion, bow_contacts })
    }

    // Reject unknown member IDs before loading material packs or reducing bodies.
    pub fn validate_members(&self, strings: usize) -> Result<()> {
        if self.entries.iter().any(|e| e.location.string >= strings
            || e.receiver_bridge.is_some_and(|i| i >= strings)) {
            return Err("stop or receiver footprint names an absent material string".into());
        }
        Ok(())
    }

    /// Reuse the runtime's complete trajectory admission before pack I/O.
    pub fn validate_performance(&self, spec: &input::Spec, gate: &CancelGate) -> Result<()> {
        poll(gate)?;
        self.validate_members(spec.strings.len())?;
        if spec.bow_ribbon.is_some() && !self.bow_contacts.is_empty() {
            return Err("uniform bow ribbon requires prescribed total load, not a bow-contact normal-load binding".into());
        }
        for binding in &self.bow_contacts {
            let track = spec.schedule.tracks().get(binding.bow)
                .ok_or("bow-contact names an absent scheduled bow")?;
            let stop = &self.entries[binding.stop];
            if fs_scenario::gesture::path::legacy_member(&track.target, "bow") != Some(stop.location.string) {
                return Err("bow-contact and stop must address the same material string".into());
            }
            for value in std::iter::once(&track.initial).chain(track.events.iter().map(|e| &e.value)) {
                if !matches!(value, GestureValue::Bow {normal_force_n, station, ..}
                    if *normal_force_n == 0.0 && station.to_bits() == stop.location.station_fraction.to_bits()) {
                    return Err("contact-loaded bow requires the stop station and no separately prescribed normal force".into());
                }
            }
        }
        if let Some(motion) = &self.motion {
            motion.validate(self.entries.len(), spec.samples, spec.rate, gate)?;
        }
        if let Some(motion) = &self.station_motion {
            let initial: Vec<_> = self.entries.iter().map(|e| e.location.station_fraction).collect();
            let lengths: Vec<_> = self.entries.iter().map(|e| spec.strings[e.location.string].length).collect();
            motion.validate(&initial, &lengths, spec.samples, spec.rate, gate)?;
            for binding in &self.bow_contacts {
                if motion.knots.iter().any(|k| k.station_fractions[binding.stop].to_bits() != initial[binding.stop].to_bits()) {
                    return Err("contact-loaded bow normals must retain their original stations".into());
                }
            }
        }
        Ok(())
    }

    pub fn motion_program(&self) -> Option<&StopMotionProgram> { self.motion.as_ref() }
    pub fn station_program(&self) -> Option<&StopStationProgram> { self.station_motion.as_ref() }
    pub fn bow_contacts(&self) -> &[BowContactLoad] { &self.bow_contacts }

    /// Bind on the already assembled and optionally preloaded instrument.
    pub fn bind_motion(&self, performance: MultiBowPerformance, gate: &CancelGate) -> Result<MultiBowPerformance> {
        let performance = match &self.motion {
            Some(motion) => performance.with_stop_motion(motion.clone(), gate)?,
            None => performance,
        };
        Ok(match &self.station_motion {
            Some(motion) => performance.with_stop_station_motion(motion.clone(), gate)?,
            None => performance,
        })
    }

    pub fn bind(&self, bank: MultiBowStringBank, gate: &CancelGate) -> Result<MultiBowStringBank> {
        poll(gate)?;
        let strings = bank.source().bank().strings();
        self.validate_members(strings.len())?;
        let mut stops = Vec::with_capacity(self.entries.len());
        for entry in &self.entries {
            poll(gate)?;
            let body = entry.receiver_bridge.map(|i| -> Result<_> {
                let links = &strings[i].connections;
                if links.len() != 1 { return Err("receiver bridge must have exactly one actual body footprint".into()); }
                let mut body = links[0].body.clone();
                if let StringStopSide::Along(normal) = entry.side {
                    // A scalar bridge footprint measures only its declared
                    // axis. Do not pretend that it measures another direction.
                    // This file format uses a common transverse frame; rotated
                    // receiver frames need explicit geometry outside this adapter.
                    let axis = match &strings[i].source {
                        MaterialStringSource::WithPolarizations(p) => p.directions().bridge(),
                        _ => [1.0, 0.0],
                    };
                    let n = normal.components();
                    let cross = n[0]*axis[1] - n[1]*axis[0];
                    let dot = n[0]*axis[0] + n[1]*axis[1];
                    if cross.abs() > 1e-12 || (dot.abs()-1.0).abs() > 1e-12 {
                        return Err("explicit stop normal must align with the selected receiver bridge axis".into());
                    }
                    for value in &mut body.shapes { *value *= dot.signum(); }
                }
                Ok(body)
            }).transpose()?;
            stops.push(MaterialStringStop { location: entry.location, side: entry.side, body,
                law: entry.law.clone(), limits: entry.limits });
        }
        // This admission precedes scheduled preload; the stops cannot disappear
        // from equilibrium and reappear only after a pluck has been initialized.
        let bank = bank.with_string_stops(stops, FrictionObstacleConfig {
            max_contacts: MAX_STOPS, max_setup_terms: 1_048_576,
        }, gate)?;
        // Actual material geometry owns co-location, perpendicular directions and
        // torsional handedness; this adapter never reconstructs its modal maps.
        Ok(bank.with_bow_contact_loads(self.bow_contacts.clone(), gate)?)
    }
}

/// Recognize each optional profile once, in either order. Paths remain OS-native.
pub fn profile_paths(args: &[std::ffi::OsString]) -> Result<(Option<&Path>, Option<&Path>)> {
    if args.len() % 2 != 0 { return Err("profile option needs a path".into()); }
    let (mut body, mut stops) = (None, None);
    for pair in args.chunks_exact(2) {
        let slot = if pair[0] == "--body-profile" { &mut body }
            else if pair[0] == "--stop-profile" { &mut stops }
            else { return Err("unknown material-string option".into()); };
        if slot.is_some() { return Err("a profile option cannot be repeated".into()); }
        *slot = Some(Path::new(&pair[1]));
    }
    Ok((body, stops))
}

// Terminal motion blocks retain the exact authored knots in stop declaration
// order. The shared runtime validates timing, units and speeds on the real clock.
// There is no alternate interpolation, contact law or source-gap mutation here.
fn parse_motion(fields: &[&str], lines: &mut std::str::Lines<'_>, stops: usize) -> Result<StopMotionProgram> {
    if fields.len() != 4 || stops == 0 {
        return Err("motion requires existing stops and MAX_OFFSET_M MAX_SPEED_M_S KNOT_COUNT".into());
    }
    let limits = ObstacleMotionLimits {
        maximum_abs_offset_m: scalar(fields[1])?, maximum_speed_m_s: scalar(fields[2])?,
    };
    limits.validate()?;
    let count: usize = fields[3].parse()?;
    if count == 0 || count > MAX_MOTION_KNOTS
        || stops.checked_mul(count).is_none_or(|n| n > MAX_MOTION_COORDINATES) {
        return Err("motion exceeds its nonempty knot/coordinate budget".into());
    }
    let mut knots = Vec::with_capacity(count);
    for _ in 0..count {
        let row: Vec<_> = lines.next().ok_or("missing motion knot")?.split_ascii_whitespace().collect();
        if row.len() != stops + 2 || row[0] != "knot" {
            return Err("expected knot ELAPSED_SAMPLE OFFSET_M for every declared stop".into());
        }
        knots.push(StopMotionKnot { sample: row[1].parse()?,
            offsets_m: row[2..].iter().map(|x| scalar(x)).collect::<Result<Vec<_>>>()? });
    }
    Ok(StopMotionProgram { knots, limits, max_knots: MAX_MOTION_KNOTS,
        max_coordinate_values: MAX_MOTION_COORDINATES })
}

fn parse_station_motion(fields: &[&str], lines: &mut std::str::Lines<'_>, stops: usize) -> Result<StopStationProgram> {
    if fields.len() != 5 || stops == 0 {
        return Err("station-motion requires existing stops and MAX_SPEED_M_S MAX_STEP_M MAX_UPDATE_TERMS KNOT_COUNT".into());
    }
    let limits = StopStationMotionLimits { maximum_speed_m_s: scalar(fields[1])?,
        maximum_step_m: scalar(fields[2])?, max_update_terms: fields[3].parse()? };
    limits.validate()?;
    if limits.max_update_terms > 1_048_576 { return Err("station geometry update work exceeds 1048576 terms".into()); }
    let count: usize = fields[4].parse()?;
    if count == 0 || count > MAX_MOTION_KNOTS
        || stops.checked_mul(count).is_none_or(|n| n > MAX_MOTION_COORDINATES) {
        return Err("station motion exceeds its nonempty knot/coordinate budget".into());
    }
    let mut knots = Vec::with_capacity(count);
    for _ in 0..count {
        let row: Vec<_> = lines.next().ok_or("missing station knot")?.split_ascii_whitespace().collect();
        if row.len() != stops + 2 || row[0] != "station-knot" {
            return Err("expected station-knot ELAPSED_SAMPLE STATION_FRACTION for every declared stop".into());
        }
        knots.push(StopStationKnot { sample: row[1].parse()?,
            station_fractions: row[2..].iter().map(|x| scalar(x)).collect::<Result<Vec<_>>>()? });
    }
    Ok(StopStationProgram { knots, limits, max_knots: MAX_MOTION_KNOTS,
        max_coordinate_values: MAX_MOTION_COORDINATES })
}

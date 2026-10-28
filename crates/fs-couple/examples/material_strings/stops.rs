//! File-driven compliant stops and normal trajectories in the original material/bow/contact performer.
//! This profile selects authored laws or explicit bulk/interface/model sources;
//! no tissue or fret properties are inferred. Existing source-derived centreline and receiver shapes own force
//! projection; the common contact, preload and time solvers own all mechanics.

use super::*;
#[path = "stops/material.rs"]
mod material;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::{
    MultiBowStringBank, stops::{MaterialStringStop, StringStopSide},
    surfaces::BowBodySurface,
    ribbon::UniformBowRibbon,
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

/// A cell reference follows the original ribbon's exact geometry projection.
#[derive(Clone, Copy, Debug)]
enum Station {
    Fraction(f64),
    RibbonCell(usize),
}
#[derive(Clone, Copy, Debug)]
enum Receiver {
    Fixed,
    Bridge(usize),
    Fiber(usize),
}
#[derive(Clone, Debug)]
enum Law {
    Authored(Obstacle),
    Material(material::Declaration),
}
#[derive(Clone, Debug)]
struct Entry {
    name: String,
    string: usize,
    station: Station,
    side: StringStopSide,
    receiver: Receiver,
    law: Law,
    limits: ModalContactConfig,
}
impl Entry {
    fn location(&self, ribbon: Option<&UniformBowRibbon>) -> Result<StringBowLocation> {
        let station_fraction = match self.station {
            Station::Fraction(value) => value,
            Station::RibbonCell(cell) => {
                let ribbon = ribbon.ok_or("ribbon:CELL requires a finite-width bow")?;
                if ribbon.centre().string != self.string {
                    return Err("ribbon cell and stop must name the same material string".into());
                }
                *ribbon.stations().get(cell).ok_or("stop names an absent ribbon cell")?
            },
        };
        Ok(StringBowLocation { string: self.string, station_fraction })
    }
}
#[derive(Clone, Debug)]
struct BowSurface {
    bow: usize,
    receiver_bridge: usize,
}
#[derive(Clone, Debug)]
pub struct Profile {
    pub bytes: Vec<u8>,
    entries: Vec<Entry>,
    motion: Option<StopMotionProgram>,
    station_motion: Option<StopStationProgram>,
    bow_contacts: Vec<BowContactLoad>,
    bow_surfaces: Vec<BowSurface>,
    sources: material::Sources,
}

fn scalar(word: &str) -> Result<f64> {
    let value: f64 = word.parse()?;
    if !value.is_finite() { return Err("stop scalar must be finite".into()); }
    Ok(value)
}
// A scalar physical footprint cannot stand in for an oblique displacement.
// Both directions use the profile's explicitly shared transverse frame.
fn surface_sign(tangent: [f64; 2], axis: [f64; 2]) -> Result<f64> {
    let cross = tangent[0]*axis[1] - tangent[1]*axis[0];
    let dot = tangent[0]*axis[0] + tangent[1]*axis[1];
    if cross.abs() > 1e-12 || (dot.abs()-1.0).abs() > 1e-12 {
        return Err("bow tangent must align with the selected receiver bridge axis".into());
    }
    Ok(dot.signum())
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
        let mut bow_surfaces: Vec<BowSurface> = Vec::new();
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
                return Err("stop, bow-contact and bow-surface declarations must precede the geometry motion blocks".into());
            }
            if fields.first() == Some(&"bow-surface") {
                if fields.len() != 3 { return Err("expected bow-surface BOW_INDEX bridge:STRING".into()); }
                let bow: usize = fields[1].parse()?;
                let receiver_bridge = fields[2].strip_prefix("bridge:")
                    .ok_or("bow surface requires an actual bridge:STRING receiver")?.parse::<usize>()?;
                if bow >= 16 || bow_surfaces.iter().any(|surface| surface.bow == bow) {
                    return Err("bow surfaces require distinct actual bows below the 16-contact bound".into());
                }
                bow_surfaces.push(BowSurface { bow, receiver_bridge });
                continue;
            }
            if fields.first() == Some(&"bow-contact") {
                if fields.len() != 3 { return Err("expected bow-contact BOW_INDEX STOP_ID".into()); }
                let bow: usize = fields[1].parse()?;
                let stop = entries.iter().position(|e| e.name == fields[2])
                    .ok_or("bow-contact must name a preceding stop declaration")?;
                if bow >= 63 || bow_contacts.iter().any(|b| b.bow == bow || b.stop == stop) {
                    return Err("bow contacts require distinct actual bows and distinct normal contacts".into());
                }
                bow_contacts.push(BowContactLoad { bow, stop });
                continue;
            }
            if entries.len() == MAX_STOPS { return Err("stop profile exceeds 32 contacts".into()); }
            if fields.len() < 5 || !matches!(fields[0], "stop" | "material-stop") {
                return Err("expected a complete stop or material-stop declaration".into());
            }
            let material = fields[0] == "material-stop";
            let width = if material { 20 } else { 12 };
            let name = fields[1];
            if name.len() > 64 || !name.bytes().all(|b| b.is_ascii_alphanumeric() || b"_-".contains(&b))
                || entries.iter().any(|e| e.name == name) {
                return Err("stop identities must be unique short ASCII tokens".into());
            }
            let string = fields[2].parse()?;
            let station = if let Some(cell) = fields[3].strip_prefix("ribbon:") {
                Station::RibbonCell(cell.parse()?)
            } else {
                Station::Fraction(scalar(fields[3])?)
            };
            let (side, values) = match fields[4] {
                "positive" if fields.len() == width => (StringStopSide::Positive, &fields[5..]),
                "negative" if fields.len() == width => (StringStopSide::Negative, &fields[5..]),
                "normal" if fields.len() == width + 2 => (StringStopSide::Along(
                    TransverseDirection::try_new([scalar(fields[5])?, scalar(fields[6])?])?), &fields[7..]),
                _ => return Err("stop side must be positive, negative or normal NX NY with a complete record".into()),
            };
            let receiver = if values[0] == "fixed" { Receiver::Fixed }
                else if let Some(index) = values[0].strip_prefix("bridge:") {
                    Receiver::Bridge(index.parse()?)
                } else if let Some(index) = values[0].strip_prefix("fiber:") {
                    Receiver::Fiber(index.parse()?)
                } else {
                    return Err("receiver must be fixed, bridge:STRING or fiber:CELL".into());
                };
            let gap = scalar(values[1])?;
            let max_force = scalar(values[values.len()-2])?;
            let max_penetration = scalar(values[values.len()-1])?;
            if matches!(station, Station::Fraction(value) if !(0.0 < value && value < 1.0))
                || matches!(station, Station::RibbonCell(cell) if cell >= 63)
                || matches!(receiver, Receiver::Fiber(cell) if cell >= 63)
                || max_force <= 0.0 || max_penetration <= 0.0 {
                return Err("stop geometry and physical limits are invalid".into());
            }
            let law = if material {
                if !matches!(receiver, Receiver::Fixed) {
                    return Err("material-stop currently requires an explicit fixed prescribed receiver".into());
                }
                Law::Material(material::Declaration::parse(values, gap, max_penetration)?)
            } else {
                let stiffness = scalar(values[2])?;
                let alpha = scalar(values[3])?;
                let loss = scalar(values[4])?;
                if stiffness <= 0.0 || alpha < 1.0 || loss < 0.0 {
                    return Err("authored stop law is invalid".into());
                }
                // Unit weight adds no fictitious area; gap is surface clearance.
                Law::Authored(Obstacle::new(vec![-1.0], 1, 1, vec![gap], vec![1.0], stiffness, alpha,
                    format!("authored string-stop profile v1: {name}; not measured contact data"))?
                    .with_internal_loss(loss)?)
            };
            entries.push(Entry { name: name.into(), string, station, side, receiver, law,
                limits: ModalContactConfig { max_iterations: 96, maximum_force_n: max_force,
                    maximum_penetration_m: max_penetration, force_absolute_tolerance_n: 1e-11,
                    force_relative_tolerance: 1e-9 } });
        }
        if entries.is_empty() && bow_surfaces.is_empty() {
            return Err("stop profile must declare at least one contact or bow surface".into());
        }
        Ok(Self { bytes: bytes.to_vec(), entries, motion, station_motion, bow_contacts, bow_surfaces,
            sources: material::Sources::default() })
    }

    // Reject unknown member IDs before loading material packs or reducing bodies.
    pub fn validate_members(&self, strings: usize) -> Result<()> {
        if self.entries.iter().any(|e| e.string >= strings
            || matches!(e.receiver, Receiver::Bridge(i) if i >= strings))
            || self.bow_surfaces.iter().any(|surface| surface.receiver_bridge >= strings) {
            return Err("stop or receiver footprint names an absent material string".into());
        }
        Ok(())
    }

    /// Reuse the runtime's complete trajectory admission before pack I/O.
    pub fn validate_performance(&self, spec: &input::Spec, gate: &CancelGate) -> Result<()> {
        poll(gate)?;
        self.validate_members(spec.strings.len())?;
        if spec.bow_ribbon.is_some() && !self.bow_surfaces.is_empty() {
            return Err("ribbon tangential surfaces are supplied by bow_fibers, not bridge footprints".into());
        }
        // The runtime's original projection owns both the centre-first cell
        // order and every station bit; input files can name cells directly.
        let ribbon = match spec.bow_ribbon {
            Some((width, cells)) => {
                let track = spec.schedule.tracks().first()
                    .ok_or("ribbon needs an original bow gesture track")?;
                let string = fs_scenario::gesture::path::legacy_member(&track.target, "bow")
                    .ok_or("ribbon needs an actual material-string bow")?;
                let GestureValue::Bow { station, .. } = track.initial else {
                    return Err("ribbon needs a complete original bow gesture".into());
                };
                let member = spec.strings.get(string)
                    .ok_or("ribbon names an absent material string")?;
                Some(UniformBowRibbon::try_new(StringBowLocation { string, station_fraction: station },
                    member.length, width, cells)?)
            },
            None => None,
        };
        let locations = self.entries.iter().map(|entry| entry.location(ribbon.as_ref()))
            .collect::<Result<Vec<_>>>()?;
        let mut fiber_normal = None;
        for (index, entry) in self.entries.iter().enumerate() {
            if let Receiver::Fiber(cell) = entry.receiver {
                let ribbon = ribbon.as_ref().ok_or("fiber receiver requires bow_ribbon")?;
                if spec.bow_fibers.is_none() || cell >= ribbon.stations().len()
                    || locations[index].string != ribbon.centre().string
                    || locations[index].station_fraction.to_bits() != ribbon.stations()[cell].to_bits()
                    || !self.bow_contacts.iter().any(|binding| binding.bow == cell && binding.stop == index) {
                    return Err("fiber receiver must match its original ribbon cell and bow-contact binding".into());
                }
                let StringStopSide::Along(normal) = entry.side else {
                    return Err("fiber contact requires an explicit normal in the shared transverse frame".into());
                };
                if fiber_normal.is_some_and(|prior| prior != normal) {
                    return Err("parallel ribbon fibers must share one contact-normal direction".into());
                }
                fiber_normal = Some(normal);
            }
        }
        for surface in &self.bow_surfaces {
            let track = spec.schedule.tracks().get(surface.bow)
                .ok_or("bow-surface names an absent scheduled bow")?;
            let string = fs_scenario::gesture::path::legacy_member(&track.target, "bow")
                .and_then(|index| spec.strings.get(index))
                .ok_or("bow-surface requires an actual material-string bow")?;
            let tangent = string.polarizations.map_or([1.0, 0.0], |directions| directions.bow());
            let axis = spec.strings[surface.receiver_bridge].polarizations
                .map_or([1.0, 0.0], |directions| directions.bridge());
            surface_sign(tangent, axis)?;
        }
        for binding in &self.bow_contacts {
            let stop = locations[binding.stop];
            let (track_index, station) = if let Some(ribbon) = &ribbon {
                let cell = ribbon.stations().get(binding.bow)
                    .ok_or("bow-contact names an absent ribbon cell")?;
                if stop.station_fraction.to_bits() != cell.to_bits() {
                    return Err("bow-contact and stop must occupy the same ribbon cell".into());
                }
                (0, ribbon.centre().station_fraction)
            } else { (binding.bow, stop.station_fraction) };
            let track = spec.schedule.tracks().get(track_index)
                .ok_or("bow-contact names an absent scheduled bow")?;
            if fs_scenario::gesture::path::legacy_member(&track.target, "bow") != Some(stop.string) {
                return Err("bow-contact and stop must address the same material string".into());
            }
            for value in std::iter::once(&track.initial).chain(track.events.iter().map(|e| &e.value)) {
                if !matches!(value, GestureValue::Bow {normal_force_n, station: value_station, ..}
                    if *normal_force_n == 0.0 && value_station.to_bits() == station.to_bits()) {
                    return Err("contact-loaded bow requires its original gesture station and zero prescribed normal force".into());
                }
            }
        }
        if let Some(motion) = &self.motion {
            motion.validate(self.entries.len(), spec.samples, spec.rate, gate)?;
            for (index, entry) in self.entries.iter().enumerate() {
                if matches!(entry.receiver, Receiver::Fiber(_))
                    && motion.knots.iter().any(|knot| knot.offsets_m[index] != 0.0) {
                    return Err("moving fiber supports require their frame-acceleration forces; use frozen fiber clearance".into());
                }
            }
        }
        if let Some(motion) = &self.station_motion {
            let initial: Vec<_> = locations.iter().map(|location| location.station_fraction).collect();
            let lengths: Vec<_> = locations.iter().map(|location| spec.strings[location.string].length).collect();
            motion.validate(&initial, &lengths, spec.samples, spec.rate, gate)?;
            for binding in &self.bow_contacts {
                if motion.knots.iter().any(|k| k.station_fractions[binding.stop].to_bits() != initial[binding.stop].to_bits()) {
                    return Err("contact-loaded normals follow the bow program; independent stop paths retain their original station entries".into());
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

    pub fn bind(&self, bank: MultiBowStringBank, spec: &input::Spec, catalog: &mut Catalog,
        loader: &mut impl FnMut(&str) -> Result<NormalizedPack>, gate: &CancelGate) -> Result<MultiBowStringBank> {
        poll(gate)?;
        // Resolve all card/model/geometry associations before binding any runtime
        // contact. A sourced record has no independently authored K/alpha/loss.
        let sources = self.sources.resolve(&self.entries, spec, catalog, loader, gate)?;
        let strings = bank.source().bank().strings();
        self.validate_members(strings.len())?;
        if bank.source().ribbon().is_some() && !self.bow_surfaces.is_empty() {
            return Err("uniform bow ribbon does not declare individual tangential receiver footprints".into());
        }
        let mut surfaces = Vec::with_capacity(self.bow_surfaces.len());
        for surface in &self.bow_surfaces {
            poll(gate)?;
            let bow = bank.source().bows().get(surface.bow)
                .ok_or("bow-surface names an absent material bow")?;
            let tangent = match &strings[bow.location.string].source {
                MaterialStringSource::WithPolarizations(p) => p.directions().bow(),
                _ => [1.0, 0.0],
            };
            let receiver = &strings[surface.receiver_bridge];
            if receiver.connections.len() != 1 {
                return Err("receiver bridge must have exactly one actual body footprint".into());
            }
            let axis = match &receiver.source {
                MaterialStringSource::WithPolarizations(p) => p.directions().bridge(),
                _ => [1.0, 0.0],
            };
            let sign = surface_sign(tangent, axis)?;
            let mut body = receiver.connections[0].body.clone();
            for value in &mut body.shapes { *value *= sign; }
            surfaces.push(BowBodySurface { bow: surface.bow, body });
        }
        let mut stops = Vec::with_capacity(self.entries.len());
        for (index, entry) in self.entries.iter().enumerate() {
            poll(gate)?;
            let body = match entry.receiver {
                Receiver::Fixed => None,
                Receiver::Bridge(i) => {
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
                Some(body)
                },
                Receiver::Fiber(cell) => {
                    let fibers = bank.source().fiber_ribbon()
                        .ok_or("fiber receiver requires retained material-fiber components")?;
                    let location = entry.location(bank.source().ribbon())?;
                    let bow = bank.source().bows().get(cell).ok_or("fiber receiver names an absent ribbon cell")?;
                    if location != bow.location || !matches!(entry.side, StringStopSide::Along(_)) {
                        return Err("fiber normal contact must occupy its original ribbon cell".into());
                    }
                    Some(fibers.normal_stop_receivers().get(cell)
                        .ok_or("fiber receiver names an absent normal map")?.clone())
                },
            };
            let law = match &entry.law {
                Law::Authored(law) => law.clone(),
                Law::Material(declaration) => {
                    let source = &sources.iter().find(|(i, _)| *i == index)
                        .ok_or("material stop source was not resolved")?.1;
                    fs_couple::unilateral_contact::binding::bound_point_obstacle(
                        source, vec![-1.0], declaration.gap_m, &entry.name)?
                },
            };
            stops.push(MaterialStringStop { location: entry.location(bank.source().ribbon())?, side: entry.side, body,
                law, limits: entry.limits });
        }
        // This admission precedes scheduled preload; the stops cannot disappear
        // from equilibrium and reappear only after a pluck has been initialized.
        let bank = bank.with_bow_body_surfaces(surfaces, gate)?;
        let bank = bank.with_string_stops(stops, FrictionObstacleConfig {
            max_contacts: MAX_STOPS, max_setup_terms: 1_048_576,
        }, gate)?.with_string_stop_sources(sources, gate)?;
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

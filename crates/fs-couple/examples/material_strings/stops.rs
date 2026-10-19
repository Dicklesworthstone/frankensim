//! File-driven fixed stops in the original material/bow/contact performer.
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

pub const MAX_BYTES: usize = 16_384;
const MAX_STOPS: usize = 32;

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
        for line in lines {
            if entries.len() == MAX_STOPS { return Err("stop profile exceeds 32 contacts".into()); }
            let fields: Vec<_> = line.split_ascii_whitespace().collect();
            if fields.len() != 12 || fields[0] != "stop" {
                return Err("expected stop ID STRING STATION SIDE RECEIVER GAP K ALPHA LOSS MAX_FORCE MAX_PENETRATION".into());
            }
            let name = fields[1];
            if name.len() > 64 || !name.bytes().all(|b| b.is_ascii_alphanumeric() || b"_-".contains(&b))
                || entries.iter().any(|e| e.name == name) {
                return Err("stop identities must be unique short ASCII tokens".into());
            }
            let location = StringBowLocation { string: fields[2].parse()?, station_fraction: scalar(fields[3])? };
            let side = match fields[4] {
                "positive" => StringStopSide::Positive,
                "negative" => StringStopSide::Negative,
                _ => return Err("stop side must be positive or negative".into()),
            };
            let receiver_bridge = if fields[5] == "fixed" { None } else {
                Some(fields[5].strip_prefix("bridge:").ok_or("receiver must be fixed or bridge:STRING")?.parse::<usize>()?)
            };
            let gap = scalar(fields[6])?;
            let stiffness = scalar(fields[7])?;
            let alpha = scalar(fields[8])?;
            let loss = scalar(fields[9])?;
            let max_force = scalar(fields[10])?;
            let max_penetration = scalar(fields[11])?;
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
        Ok(Self { bytes: bytes.to_vec(), entries })
    }

    // Reject unknown member IDs before loading material packs or reducing bodies.
    pub fn validate_members(&self, strings: usize) -> Result<()> {
        if self.entries.iter().any(|e| e.location.string >= strings
            || e.receiver_bridge.is_some_and(|i| i >= strings)) {
            return Err("stop or receiver footprint names an absent material string".into());
        }
        Ok(())
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
                Ok(links[0].body.clone())
            }).transpose()?;
            stops.push(MaterialStringStop { location: entry.location, side: entry.side, body,
                law: entry.law.clone(), limits: entry.limits });
        }
        // This admission precedes scheduled preload; the stops cannot disappear
        // from equilibrium and reappear only after a pluck has been initialized.
        Ok(bank.with_string_stops(stops, FrictionObstacleConfig {
            max_contacts: MAX_STOPS, max_setup_terms: 1_048_576,
        }, gate)?)
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

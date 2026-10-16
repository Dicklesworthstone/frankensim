//! Bounded geometry/input adapter; material laws and gesture interpolation live elsewhere.
use super::*;
use std::str::{FromStr, Lines, SplitAsciiWhitespace};
use fs_material::gas_species::input::read_gas_state;
use fs_scenario::gesture::GestureValue;

pub const MAX_INPUT_BYTES: usize = 1_048_576;
pub const MAX_PACK_BYTES: usize = 8 * 1_048_576;
pub const MAX_TOTAL_PACK_BYTES: usize = 64 * 1_048_576;
pub const MAX_MODES: usize = 128;
pub const MAX_EVENTS: usize = 16_384;
pub const MAX_CONTROLS: usize = 262_144;
pub const MAX_WORK: u64 = 16_777_216;

#[derive(Clone, Debug)]
pub struct Wire {
    pub pack: String, pub length: f64, pub radius: f64, pub prestress: StringPrestress,
    pub section: WireSection,
    pub bending_modes: usize, pub twist_modes: usize, pub bridge: SharedPlateBridge,
}
/// V2 describes perfectly bonded circular layers, never helical winding/slip.
#[derive(Clone, Debug)]
pub enum WireSection {
    Homogeneous,
    Bonded { core_radius_m: f64, sleeves: Vec<WireLayer> },
}
#[derive(Clone, Debug)]
pub struct WireLayer { pub pack: String, pub outer_radius_m: f64 }
const MAX_LAYERS: usize = 16;
#[derive(Clone, Debug)]
pub struct Body {
    pub pack: String, pub length: f64, pub width: f64, pub thickness: f64,
    pub modes: usize, pub listener: f64,
}
#[derive(Clone, Debug)]
pub struct Spec {
    pub bytes: Vec<u8>, pub rate: u32, pub samples: u64, pub full_scale: f64, pub block: usize,
    pub temperature: f64, pub omega: f64, pub gas: GasState, pub body: Body, pub strings: Vec<Wire>,
    pub cavity: Option<MaterialBodyCavity>, pub law: FrictionLaw, pub forces: Option<StringForceProgram>,
    pub schedule: GestureSchedule,
}
impl Spec {
    pub fn parse(bytes: &[u8]) -> Result<Self> {
        if bytes.len() > MAX_INPUT_BYTES { return Err("material-string input exceeds 1 MiB".into()); }
        let text = std::str::from_utf8(bytes)?;
        let (header, schedule) = text.split_once("\nschedule\n")
            .ok_or("expected LF-delimited schedule followed by canonical gesture bytes")?;
        let mut r = Reader { lines: header.lines(), line: 0 };
        // V1 retains its original whitespace acceptance and record semantics.
        let mut version = r.raw_row()?;
        let layered = match version.word()? {
            "frankensim-material-strings-v1" => false,
            "frankensim-material-strings-v2" => true,
            _ => return Err("unsupported material-string input version".into()),
        };
        version.end()?;
        let mut row = r.row("audio")?;
        let rate: u32 = row.parse()?; let samples: u64 = row.parse()?;
        let full_scale = row.positive()?; let block = row.count(65_536)?; row.end()?;
        if rate == 0 || rate > 192_000 || samples == 0 || samples > 28_800_000
            || samples > u64::from(rate)*600 || block == 0 {
            return Err("audio requires a bounded positive clock, horizon and callback".into());
        }
        let mut row = r.row("state")?;
        let temperature = row.positive()?; let omega = row.scalar()?; row.end()?;
        if omega < 0.0 { return Err("source angular frequency must be nonnegative".into()); }
        let mut row = r.row("gas")?;
        let gas = read_gas_state(&mut row.words, None)?; row.end()?;
        let mut row = r.row("body")?;
        let body = Body { pack: row.path()?, length: row.positive()?, width: row.positive()?,
            thickness: row.positive()?, modes: row.count(16)?, listener: row.positive()? };
        row.end()?;
        let mut row = r.row("strings")?; let count = row.count(8)?; row.end()?;
        if count == 0 || body.modes == 0 { return Err("nonempty material strings and body modes required".into()); }
        let mut strings = Vec::with_capacity(count); let mut modes = body.modes;
        for _ in 0..count {
            let mut row = r.row("string")?;
            let pack = row.path()?; let length = row.positive()?; let radius = row.positive()?;
            let prestress = if layered {
                match row.word()? {
                    "tension" => StringPrestress::FixedTension(row.positive()?),
                    "extension" => StringPrestress::FixedExtension {
                        stress_free_length_m: row.positive()?, linear_strain_limit: row.positive()?,
                    },
                    _ => return Err("v2 support must be tension N or extension L0 STRAIN_LIMIT".into()),
                }
            } else { StringPrestress::FixedTension(row.positive()?) };
            let mut wire = Wire { pack, length, radius, prestress, section: WireSection::Homogeneous,
                bending_modes: row.count(MAX_MODES)?, twist_modes: row.count(MAX_MODES)?,
                bridge: SharedPlateBridge { station_fraction: row.scalar()?, stiffness_n_m: row.scalar()?,
                    damping_n_s_m: row.scalar()?, rest_extension_m: 0.0 } };
            row.end()?;
            if layered { wire.section = r.section(radius)?; }
            if wire.bending_modes == 0 || !(0.0 < wire.bridge.station_fraction && wire.bridge.station_fraction < 1.0)
                || wire.bridge.stiffness_n_m < 0.0 || wire.bridge.damping_n_s_m < 0.0
                || !(wire.bridge.stiffness_n_m > 0.0 || wire.bridge.damping_n_s_m > 0.0) {
                return Err("string needs transverse modes and a physical interior bridge connection".into());
            }
            modes += wire.bending_modes + wire.twist_modes; strings.push(wire);
        }
        let mut row = r.row("cavity")?;
        let cavity = match row.word()? {
            "none" => None,
            "neck" => {
                modes += 1;
                let description = HelmholtzCavity { volume_m3: row.positive()?, neck_radius_m: row.positive()?,
                    neck_length_m: row.positive()? };
                let resistance = row.scalar()?; let listener_m = row.positive()?;
                Some(MaterialBodyCavity { cavity: HelmholtzModalCavity::new(description, gas,
                    CavityNeckLoss::Resistance(resistance))?, modal_budget: modal_budget(), listener_m })
            },
            _ => return Err("cavity must be none or neck V R L RESISTANCE LISTENER_M".into()),
        };
        row.end()?;
        if modes > MAX_MODES { return Err("complete instrument exceeds 128 modes".into()); }
        let mut row = r.row("friction")?;
        let law = match row.word()? {
            "coulomb" => FrictionLaw::Coulomb { static_mu: row.scalar()?, kinetic_mu: row.scalar()? },
            "stribeck" => FrictionLaw::Stribeck { static_mu: row.scalar()?, kinetic_mu: row.scalar()?,
                characteristic_speed: row.positive()?, viscous_per_speed: row.scalar()? },
            _ => return Err("unsupported physical friction law".into()),
        };
        row.end()?; law.regularized_traction_1d(0.0, 0.0, 1.0)?;
        let mut row = r.row("forces")?;
        let initialization = match row.word()? {
            "retain" => ForceInitialization::RetainState,
            "preload" => ForceInitialization::StaticPreload,
            _ => return Err("force initialization must be retain or preload".into()),
        };
        let ports = row.count(64)?; let event_count = row.count(MAX_EVENTS)?; row.end()?;
        let mut locations = Vec::with_capacity(ports); let mut initial_force_n = Vec::with_capacity(ports);
        for _ in 0..ports {
            let mut row = r.row("port")?;
            let string = row.count(count-1)?; let station_fraction = row.scalar()?;
            let coupling = match row.word()? {
                "centreline" => StringForceCoupling::Centreline,
                "surface" => StringForceCoupling::Surface,
                _ => return Err("force coupling must be centreline or surface".into()),
            };
            if !(0.0 < station_fraction && station_fraction < 1.0) { return Err("load station must be interior".into()); }
            locations.push(StringForcePort { location: StringBowLocation { string, station_fraction }, coupling });
            initial_force_n.push(row.scalar()?); row.end()?;
        }
        let mut events = Vec::with_capacity(event_count);
        for _ in 0..event_count {
            let mut row = r.row("force")?;
            let sample: u64 = row.parse()?; let port = row.count(64)?; let force_n = row.scalar()?; row.end()?;
            if sample >= samples || port >= ports { return Err("force event exceeds its physical ports or observed window".into()); }
            events.push(ModalForceEvent { sample, voice: 0, port, force_n });
        }
        if r.lines.next().is_some() { return Err("unexpected header record".into()); }
        if ports == 0 && (event_count != 0 || initialization == ForceInitialization::StaticPreload) {
            return Err("preload and force events require an explicit load port".into());
        }
        let forces = (ports != 0).then_some(StringForceProgram { ports: locations, initial_force_n, events,
            initialization, limits: ForceRenderConfig { sample_rate_hz: rate, max_block: block,
                max_events: MAX_EVENTS, max_controls: MAX_CONTROLS, max_projection_terms: 16_777_216 } });
        // The canonical decoder reserves from declared counts. Bound them before
        // invoking it, then require exact re-encoding so nothing is silently lost.
        let mut total = 0usize; let lines = schedule.lines().count();
        for line in schedule.lines() {
            if let Some(n) = line.strip_prefix("tracks\t") {
                let n: usize = n.parse()?;
                if n == 0 || n > 16 || n > lines { return Err("bow track count exceeds source/16-port bound".into()); }
            }
            if let Some(n) = line.strip_prefix("events\t") {
                let n: usize = n.parse()?;
                if n > MAX_EVENTS || n > lines { return Err("gesture event count exceeds source bound".into()); }
                total = total.checked_add(n).ok_or("gesture count overflow")?;
                if total > MAX_EVENTS { return Err("total gesture count exceeds 16384".into()); }
            }
        }
        let schedule = GestureSchedule::from_canonical_bytes(schedule.as_bytes())?;
        if schedule.to_canonical_bytes() != text.split_once("\nschedule\n").unwrap().1.as_bytes()
            || schedule.control_rate_hz > rate || schedule.tracks().is_empty() || schedule.tracks().len() > 16 {
            return Err("canonical bow schedule and control rate no greater than physics required".into());
        }
        let last_tick = (samples-1)*u64::from(schedule.control_rate_hz)/u64::from(rate);
        let last_time = last_tick as f64/f64::from(schedule.control_rate_hz);
        for track in schedule.tracks() {
            let string = fs_scenario::gesture::path::legacy_member(&track.target,"bow").ok_or("only physical bow tracks are supported")?;
            if string >= strings.len() { return Err("bow track names an absent material string".into()); }
            let GestureValue::Bow { station, .. } = track.initial else { return Err("complete initial bow value required".into()); };
            for event in &track.events {
                if event.time_s > last_time { return Err("bow command starts after the last observed control tick".into()); }
                if !matches!(event.value, GestureValue::Bow { station: s, .. } if s.to_bits() == station.to_bits()) {
                    return Err("bow stations remain fixed; events may change speed and load, not geometry".into());
                }
            }
        }
        Ok(Self { bytes: bytes.to_vec(), rate, samples, full_scale, block, temperature, omega,
            gas, body, strings, cavity, law, forces, schedule })
    }
}

struct Reader<'a> { lines: Lines<'a>, line: usize }
struct Row<'a> { words: SplitAsciiWhitespace<'a>, line: usize }
impl<'a> Reader<'a> {
    fn raw_row(&mut self) -> Result<Row<'a>> {
        self.line += 1;
        let words = self.lines.next().ok_or_else(|| format!("missing record on line {}", self.line))?
            .split_ascii_whitespace();
        Ok(Row { words, line: self.line })
    }
    fn section(&mut self, outer_radius: f64) -> Result<WireSection> {
        let mut row = self.row("section")?;
        match row.word()? {
            "homogeneous" => { row.end()?; Ok(WireSection::Homogeneous) }
            "bonded" => {
                let core_radius_m = row.positive()?;
                let count = row.count(MAX_LAYERS - 1)?; row.end()?;
                if count == 0 || core_radius_m >= outer_radius {
                    return Err("bonded section requires a smaller solid core and 1..=15 sleeves".into());
                }
                let mut sleeves = Vec::with_capacity(count);
                let mut previous = core_radius_m;
                for _ in 0..count {
                    let mut row = self.row("layer")?;
                    let pack = row.path()?; let outer_radius_m = row.positive()?; row.end()?;
                    if outer_radius_m <= previous || outer_radius_m > outer_radius {
                        return Err("bonded radii must increase strictly within the stated outer radius".into());
                    }
                    previous = outer_radius_m;
                    sleeves.push(WireLayer { pack, outer_radius_m });
                }
                if previous.to_bits() != outer_radius.to_bits() {
                    return Err("last bonded layer must end at the string's stated outer radius".into());
                }
                Ok(WireSection::Bonded { core_radius_m, sleeves })
            }
            _ => Err("section must be homogeneous or bonded CORE_RADIUS SLEEVE_COUNT".into()),
        }
    }
    fn row(&mut self, key: &str) -> Result<Row<'a>> {
        self.line += 1;
        let mut words = self.lines.next().ok_or_else(|| format!("missing {key} on line {}", self.line))?.split_ascii_whitespace();
        if words.next() != Some(key) { return Err(format!("expected {key} on line {}", self.line).into()); }
        Ok(Row { words, line: self.line })
    }
}
impl<'a> Row<'a> {
    fn word(&mut self) -> Result<&'a str> { self.words.next().ok_or_else(|| format!("missing field on line {}", self.line).into()) }
    fn parse<T: FromStr>(&mut self) -> Result<T> {
        let line = self.line;
        self.word()?.parse().map_err(|_| format!("invalid field on line {line}").into())
    }
    fn scalar(&mut self) -> Result<f64> {
        let n: f64 = self.parse()?;
        if !n.is_finite() { return Err(format!("nonfinite scalar on line {}", self.line).into()); } Ok(n)
    }
    fn positive(&mut self) -> Result<f64> {
        let n = self.scalar()?; if n <= 0.0 { return Err(format!("positive scalar required on line {}", self.line).into()); } Ok(n)
    }
    fn count(&mut self, cap: usize) -> Result<usize> {
        let n = self.parse()?; if n > cap { return Err(format!("count exceeds {cap} on line {}", self.line).into()); } Ok(n)
    }
    fn path(&mut self) -> Result<String> {
        let name = self.word()?;
        if name.len() > 512 || !std::path::Path::new(name).components()
            .all(|p| matches!(p, std::path::Component::Normal(_))) {
            return Err("material pack paths must stay relative to the performance directory".into());
        }
        Ok(name.to_string())
    }
    fn end(mut self) -> Result<()> {
        if self.words.next().is_some() { return Err(format!("extra field on line {}", self.line).into()); } Ok(())
    }
}

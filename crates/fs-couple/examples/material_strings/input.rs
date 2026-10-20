//! Bounded geometry/input adapter; material laws and gesture interpolation live elsewhere.
use super::*;
use std::str::{FromStr, Lines, SplitAsciiWhitespace};
use fs_material::gas_species::input::{GasRecord, ThermochemicalSources, read_gas_record};
use fs_scenario::gesture::GestureValue;
use fs_couple::bowed_string::runtime::reciprocal::bank::polarization::{TransverseDirection, TransverseDirections};

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
    /// With fixed extension, L0 is stress-free at this explicit reference T.
    /// None preserves the original same-temperature mechanical prescription.
    pub thermal_reference_kelvin: Option<f64>,
    pub bending_modes: usize, pub twist_modes: usize, pub bridge: SharedPlateBridge,
    pub polarizations: Option<TransverseDirections>,
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
const MAX_ROTATIONAL_BRIDGES: usize = 16;
/// A scalar rotational joint referencing actual body-profile point footprints.
#[derive(Clone, Debug, PartialEq)]
pub struct RotationalBridge {
    pub string: usize, pub station_fraction: f64,
    pub positive_bridge_string: usize, pub negative_bridge_string: usize,
    pub stiffness_n_m_per_rad: f64, pub damping_n_m_s_per_rad: f64, pub rest_angle_rad: f64,
}
#[derive(Clone, Debug)]
pub struct Body {
    pub pack: String, pub length: f64, pub width: f64, pub thickness: f64,
    pub modes: usize, pub listener: f64,
}
#[derive(Clone, Debug)]
pub struct Spec {
    pub bytes: Vec<u8>, pub rate: u32, pub samples: u64, pub full_scale: f64, pub block: usize,
    pub temperature: f64, pub omega: f64, pub gas: GasState, pub body: Body, pub strings: Vec<Wire>,
    /// Retain original thermochemical composition and pack authority with the coefficients.
    pub gas_record: GasRecord,
    pub rotational_bridges: Vec<RotationalBridge>,
    pub cylinder_gas_loading: Option<CylinderFlowResolution>,
    pub prepared_linear: Option<PreparedModalBudget>,
    pub geometric_stretching: Option<(f64, NonlinearModalConfig)>,
    pub cavity: Option<MaterialBodyCavity>, pub friction: bow_interfaces::Inputs,
    pub bow_ribbon: Option<(f64, usize)>, pub forces: Option<StringForceProgram>,
    pub schedule: GestureSchedule,
}
impl Spec {
    pub fn parse(bytes: &[u8]) -> Result<Self> {
        Self::parse_with_gas_sources(bytes, None)
    }

    pub fn parse_with_gas_sources(bytes: &[u8], sources: Option<&ThermochemicalSources<'_>>) -> Result<Self> {
        if bytes.len() > MAX_INPUT_BYTES { return Err("material-string input exceeds 1 MiB".into()); }
        let text = std::str::from_utf8(bytes)?;
        let (header, schedule) = text.split_once("\nschedule\n")
            .ok_or("expected LF-delimited schedule followed by canonical gesture bytes")?;
        let mut r = Reader { lines: header.lines(), line: 0 };
        // V1 retains its original whitespace acceptance and record semantics.
        let mut version = r.raw_row()?;
        let format = match version.word()? {
            "frankensim-material-strings-v1" => 1,
            "frankensim-material-strings-v2" => 2,
            "frankensim-material-strings-v3" => 3,
            _ => return Err("unsupported material-string input version".into()),
        };
        let layered = format >= 2;
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
        let gas_record = read_gas_record(&mut row.words, None, sources)?; row.end()?;
        let gas = *gas_record.state();
        let mut row = r.row("body")?;
        let body = Body { pack: row.path()?, length: row.positive()?, width: row.positive()?,
            thickness: row.positive()?, modes: row.count(16)?, listener: row.positive()? };
        row.end()?;
        let mut row = r.row("strings")?; let count = row.count(8)?; row.end()?;
        if count == 0 || body.modes == 0 { return Err("nonempty material strings and body modes required".into()); }
        let mut strings = Vec::with_capacity(count);
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
                thermal_reference_kelvin: None,
                polarizations: None,
                bending_modes: row.count(MAX_MODES)?, twist_modes: row.count(MAX_MODES)?,
                bridge: SharedPlateBridge { station_fraction: row.scalar()?, stiffness_n_m: row.scalar()?,
                    damping_n_s_m: row.scalar()?, rest_extension_m: 0.0 } };
            row.end()?;
            if layered { wire.section = r.section(radius)?; }
            if format >= 3 {
                let mut row = r.row("polarizations")?;
                wire.polarizations = match row.word()? {
                    "single" => None,
                    "two" => Some(TransverseDirections::try_new(
                        [row.scalar()?, row.scalar()?], [row.scalar()?, row.scalar()?])?),
                    _ => return Err("polarizations must be single or two BOW_X BOW_Y BRIDGE_X BRIDGE_Y".into()),
                };
                row.end()?;
            }
            if wire.bending_modes == 0 || !(0.0 < wire.bridge.station_fraction && wire.bridge.station_fraction < 1.0)
                || wire.bridge.stiffness_n_m < 0.0 || wire.bridge.damping_n_s_m < 0.0
                || !(wire.bridge.stiffness_n_m > 0.0 || wire.bridge.damping_n_s_m > 0.0) {
                return Err("string needs transverse modes and a physical interior bridge connection".into());
            }
            strings.push(wire);
        }
        // Optional source-backed eigenstrain at fixed supports. The original
        // extension record supplies reference L0 and the linear strain limit.
        if r.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next()) == Some("thermal-supports") {
            let mut row = r.row("thermal-supports")?;
            let supports = row.count(count)?; row.end()?;
            for _ in 0..supports {
                let mut row = r.row("thermal-support")?;
                let string = row.count(count - 1)?; let reference = row.positive()?; row.end()?;
                let wire = &mut strings[string];
                if wire.thermal_reference_kelvin.is_some()
                    || !matches!(wire.prestress, StringPrestress::FixedExtension { .. })
                    || !matches!(&wire.section, WireSection::Homogeneous) {
                    return Err("thermal support requires one declaration per homogeneous fixed-extension string".into());
                }
                wire.thermal_reference_kelvin = Some(reference);
            }
        }
        // Optional in every version. The following cavity row is left
        // untouched when absent; no new default joint or body rotation is inferred.
        let mut rotational_bridges = Vec::new();
        if r.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next()) == Some("rotational-bridges") {
            let mut row = r.row("rotational-bridges")?;
            let joints = row.count(MAX_ROTATIONAL_BRIDGES)?; row.end()?;
            rotational_bridges.reserve(joints);
            for _ in 0..joints {
                let mut row = r.row("rotational-bridge")?;
                let joint = RotationalBridge {
                    string: row.count(count - 1)?, station_fraction: row.scalar()?,
                    positive_bridge_string: row.count(count - 1)?, negative_bridge_string: row.count(count - 1)?,
                    stiffness_n_m_per_rad: row.scalar()?, damping_n_m_s_per_rad: row.scalar()?,
                    rest_angle_rad: row.scalar()?,
                };
                row.end()?;
                if strings[joint.string].twist_modes == 0
                    || !(0.0 < joint.station_fraction && joint.station_fraction < 1.0)
                    || joint.positive_bridge_string == joint.negative_bridge_string
                    || joint.stiffness_n_m_per_rad < 0.0 || joint.damping_n_m_s_per_rad < 0.0
                    || (joint.stiffness_n_m_per_rad == 0.0 && joint.damping_n_m_s_per_rad == 0.0) {
                    return Err("rotational bridge needs sourced twist, an interior station, distinct body footprints and a nonzero passive law".into());
                }
                rotational_bridges.push(joint);
            }
        }
        let cylinder_gas_loading = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("cylinder-gas-loading") {
            let mut row = r.row("cylinder-gas-loading")?;
            let resolution = CylinderFlowResolution {
                minimum_wavenumber_radius: row.positive()?,
                maximum_wavenumber_radius: row.positive()?,
                intervals: row.count(4096)?,
            };
            row.end()?;
            resolution.validate()?;
            Some(resolution)
        } else { None };
        let prepared_linear = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("prepared-linear") {
            let mut row = r.row("prepared-linear")?;
            let budget = PreparedModalBudget {
                max_modes: row.parse()?,
                max_dense_bytes: row.parse()?,
                max_multiply_terms: row.parse()?,
            };
            row.end()?;
            budget.validate()?;
            Some(budget)
        } else { None };
        let geometric_stretching = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("geometric-stretching") {
            let mut row = r.row("geometric-stretching")?;
            let slope = row.positive()?;
            let config = NonlinearModalConfig {
                max_iterations: row.count(128)?, maximum_generalized_force: row.positive()?,
                force_absolute_tolerance: row.positive()?, force_relative_tolerance: row.positive()?,
            };
            row.end()?;
            config.validate()?;
            fs_couple::render::plate::impact::string::StringStretching {
                axial_rigidity_n: 0.0, maximum_slope: slope,
            }.validate()?;
            Some((slope, config))
        } else { None };
        let mut row = r.row("cavity")?;
        let cavity = match row.word()? {
            "none" => None,
            "neck" => {
                // The `neck` record predates DEF-09: a sound hole through the
                // body's own plate, declared flanged at both ends.
                let description = HelmholtzCavity { volume_m3: row.positive()?, neck_radius_m: row.positive()?,
                    neck_length_m: row.positive()?, inner_end: fs_scenario::NeckEnd::Flanged,
                    outer_end: fs_scenario::NeckEnd::Flanged };
                let resistance = row.scalar()?; let listener_m = row.positive()?;
                Some(MaterialBodyCavity { cavity: HelmholtzModalCavity::new(description, gas,
                    CavityNeckLoss::Resistance(resistance))?, modal_budget: modal_budget(), listener_m })
            },
            "physical-neck" => {
                let volume_m3 = row.positive()?;
                let neck_radius_m = row.positive()?;
                let neck_length_m = row.positive()?;
                let mut end = || -> Result<fs_scenario::NeckEnd> {
                    match row.word()? {
                        "flanged" => Ok(fs_scenario::NeckEnd::Flanged),
                        "unflanged" => Ok(fs_scenario::NeckEnd::Unflanged),
                        _ => Err("physical neck ends must be flanged or unflanged".into()),
                    }
                };
                let inner_end = end()?;
                let outer_end = end()?;
                let listener_m = row.positive()?;
                let description = HelmholtzCavity {
                    volume_m3, neck_radius_m, neck_length_m, inner_end, outer_end,
                };
                // One coherent gas supplies inertia, compression, Bessel neck
                // loss and the declared exterior end's radiation resistance.
                // The original owner counts each end correction exactly once.
                Some(MaterialBodyCavity { cavity: HelmholtzModalCavity::new(description, gas,
                    CavityNeckLoss::RadiationAtResonance)?, modal_budget: modal_budget(), listener_m })
            },
            _ => return Err("cavity must be none, neck V R L RESISTANCE LISTENER_M, or physical-neck V R L INNER_END OUTER_END LISTENER_M".into()),
        };
        row.end()?;
        let mut row = r.row("friction")?;
        let friction = match row.word()? {
            "coulomb" => {
                let law=FrictionLaw::Coulomb { static_mu:row.scalar()?,kinetic_mu:row.scalar()? };
                row.end()?; bow_interfaces::Inputs::research(law)?
            },
            "stribeck" => {
                let law=FrictionLaw::Stribeck { static_mu:row.scalar()?,kinetic_mu:row.scalar()?,
                    characteristic_speed:row.positive()?,viscous_per_speed:row.scalar()? };
                row.end()?; bow_interfaces::Inputs::research(law)?
            },
            "interfaces" => {
                let count=row.count(16)?; row.end()?;
                if count==0 { return Err("source-bound friction requires one interface per bow track".into()); }
                let mut inputs=Vec::with_capacity(count);
                for _ in 0..count {
                    let mut row=r.row("bow-interface")?;
                    let track=row.word()?.to_string(); let pack=row.path()?;
                    let temperature=row.positive()?; let area=row.positive()?;
                    let pressure=row.positive()?; let speed=row.scalar()?;
                    let model=match row.word()? {
                        "coulomb"=>bow_interfaces::Model::Coulomb,
                        "stribeck"=>bow_interfaces::Model::Stribeck,
                        _=>return Err("bow interface model must be coulomb or stribeck".into()),
                    };
                    let names=(0..model.property_count()).map(|_|row.word().map(str::to_string))
                        .collect::<Result<Vec<_>>>()?;
                    row.end()?;
                    inputs.push(bow_interfaces::Declaration {track,pack,temperature,area,pressure,speed,model,names});
                }
                bow_interfaces::Inputs::interfaces(inputs)?
            },
            _ => return Err("unsupported physical friction law".into()),
        };
        let bow_ribbon = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("bow_ribbon") {
            let mut row = r.row("bow_ribbon")?;
            let width = row.positive()?; let cells = row.count(63)?; row.end()?;
            if cells == 0 || cells % 2 == 0 {
                return Err("bow ribbon requires an odd spatial cell count in 1..=63".into());
            }
            Some((width, cells))
        } else { None };
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
                "along" => {
                    let direction = TransverseDirection::try_new([row.scalar()?, row.scalar()?])?;
                    if strings[string].polarizations.is_none() && direction.components()[1] != 0.0 {
                        return Err("second-plane force requires polarizations two".into());
                    }
                    StringForceCoupling::CentrelineAlong(direction)
                },
                _ => return Err("force coupling must be centreline, surface or along DX DY".into()),
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
        if geometric_stretching.is_some() && initialization == ForceInitialization::StaticPreload {
            return Err("geometric stretching requires retained initial state; nonlinear static preload is not implemented".into());
        }
        if schedule.to_canonical_bytes() != text.split_once("\nschedule\n").unwrap().1.as_bytes()
            || schedule.control_rate_hz > rate || schedule.tracks().is_empty() || schedule.tracks().len() > 16 {
            return Err("canonical bow schedule and control rate no greater than physics required".into());
        }
        let last_tick = (samples-1)*u64::from(schedule.control_rate_hz)/u64::from(rate);
        let last_time = last_tick as f64/f64::from(schedule.control_rate_hz);
        if bow_ribbon.is_some() && schedule.tracks().len() != 1 {
            return Err("a uniform bow ribbon requires exactly one original total-load track".into());
        }
        for track in schedule.tracks() {
            let string = fs_scenario::gesture::path::legacy_member(&track.target,"bow").ok_or("only physical bow tracks are supported")?;
            if string >= strings.len() { return Err("bow track names an absent material string".into()); }
            let GestureValue::Bow { station, .. } = track.initial else { return Err("complete initial bow value required".into()); };
            if let Some((width, _)) = bow_ribbon {
                let half = 0.5 * width / strings[string].length;
                if !(0.0 < station - half && station - half < station
                    && station < station + half && station + half < 1.0) {
                    return Err("the complete bow ribbon must lie strictly inside its material string".into());
                }
            }
            for event in &track.events {
                if event.time_s > last_time { return Err("bow command starts after the last observed control tick".into()); }
                if !matches!(event.value, GestureValue::Bow { station: s, .. } if s.to_bits() == station.to_bits()) {
                    return Err("bow stations remain fixed; events may change speed and load, not geometry".into());
                }
            }
        }
        friction.validate_tracks(&schedule)?;
        let spec = Self { bytes: bytes.to_vec(), rate, samples, full_scale, block, temperature, omega,
            gas, gas_record, body, strings, rotational_bridges, cylinder_gas_loading, prepared_linear, geometric_stretching, cavity, friction, bow_ribbon, forces, schedule };
        spec.network_dimensions(false)?;
        Ok(spec)
    }

    pub fn bow_law(&self,port:usize)->Result<&FrictionLaw> {
        let track=self.schedule.tracks().get(port).ok_or("unknown bow contact")?;
        self.friction.law(&track.id)
    }

    pub fn bow_contact_count(&self)->usize {
        self.bow_ribbon.map_or(self.schedule.tracks().len(), |(_, cells)| cells)
    }

    /// Complete numerical allocation before source I/O. A tensor body can need
    /// one damping column per retained body mode; explicit zero viscosity may
    /// use fewer. The allowance never drops requested modes or fluid branches.
    pub fn network_dimensions(&self, tensor_body: bool) -> Result<(usize, usize, usize)> {
        self.network_dimensions_with_body_connections(tensor_body, 0)
    }

    /// Additional source-declared shell interfaces share the original complete
    /// connection/work ceilings. Reserve every declaration before pack I/O;
    /// explicit zero viscosity can reduce the exact admitted runtime count.
    pub fn network_dimensions_with_body_connections(&self, tensor_body: bool,
        additional_body_connections: usize) -> Result<(usize, usize, usize)>
    {
        let mut modes = self.body.modes.checked_add(usize::from(self.cavity.is_some()))
            .ok_or("complete mode count overflow")?;
        let mut transverse = 0usize;
        for wire in &self.strings {
            let planes = if wire.polarizations.is_some() { 2 } else { 1 };
            let bending = wire.bending_modes.checked_mul(planes)
                .ok_or("transverse mode count overflow")?;
            transverse = transverse.checked_add(bending).ok_or("transverse mode count overflow")?;
            modes = modes.checked_add(bending).and_then(|n| n.checked_add(wire.twist_modes))
                .ok_or("complete string/body mode count overflow")?;
        }
        let cavity_links = self.cavity.as_ref().map_or(0, |c|
            1 + usize::from(c.cavity.neck_resistance_n_s_m() != 0.0));
        let mut links = self.strings.len().checked_add(self.rotational_bridges.len())
            .and_then(|n| n.checked_add(cavity_links))
            .and_then(|n| n.checked_add(if tensor_body { self.body.modes } else { 0 }))
            .and_then(|n| n.checked_add(additional_body_connections))
            .ok_or("complete connection count overflow")?;
        if let Some(resolution) = self.cylinder_gas_loading {
            resolution.validate()?;
            let memory = transverse.checked_mul(resolution.intervals)
                .ok_or("fluid memory dimension overflow")?;
            modes = modes.checked_add(memory).ok_or("complete fluid mode count overflow")?;
            links = links.checked_add(memory).ok_or("complete fluid connection count overflow")?;
        }
        if let Some(budget) = self.prepared_linear {
            budget.validate()?;
            if modes > budget.max_modes {
                return Err("complete instrument exceeds the explicit prepared mode budget".into());
            }
        }
        if modes == 0 || modes > MAX_MODES {
            return Err("complete instrument, including gas memory, exceeds 128 modes".into());
        }
        if links > MAX_MODAL_CONNECTIONS {
            return Err("complete instrument exceeds the shared modal connection ceiling".into());
        }
        // The original coupled owner preflights N*(L+1)^2, not merely N*L.
        let terms = links.checked_add(1).and_then(|l| l.checked_mul(l))
            .and_then(|square| modes.checked_mul(square))
            .ok_or("complete network setup work overflow")?;
        if terms > MAX_WORK as usize {
            return Err("complete instrument exceeds the explicit network setup work budget".into());
        }
        Ok((modes, links, terms))
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

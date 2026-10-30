//! Bounded geometry/input adapter; material laws and gesture interpolation live elsewhere.
use super::*;
use std::str::{FromStr, Lines, SplitAsciiWhitespace};
use fs_material::gas_species::input::{GasRecord, ThermochemicalSources, read_gas_record};
use fs_material::relaxation::RelaxationState;
use fs_couple::bowed_string::runtime::reciprocal::bank::cavity::ShellCavityMolecularOptions;
use fs_scenario::gesture::GestureValue;
use fs_couple::bowed_string::runtime::reciprocal::bank::polarization::{TransverseDirection, TransverseDirections};
use fs_scenario::acoustic::BendingRelaxationProperties;
use fs_couple::material_assembly::measured_string::MeasuredStringProperties;
use fs_couple::bowed_string::runtime::reciprocal::bank::radiation::rayleigh::PlateRayleighOptions;
use fs_couple::bowed_string::runtime::reciprocal::bank::radiation::bem::ShellBemOptions;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::observation::ShellObserverOptions;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::observation::plate::PlateObserverOptions;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::observation::one_way::{
    OneWayRule, OneWayShellObserverOptions};

#[path = "bow_motion.rs"]
pub mod bow_motion;
#[path = "input/forces.rs"]
mod forces;

pub const MAX_INPUT_BYTES: usize = 1_048_576;
pub const MAX_PACK_BYTES: usize = 8 * 1_048_576;
pub const MAX_TOTAL_PACK_BYTES: usize = 64 * 1_048_576;
pub const MAX_MODES: usize = 128;
/// Retained body modes a `body` record may request: the shell owner's own ceiling, `MAX_SHELL_MODES`. Flat plate
/// bodies keep [`MAX_PLATE_BODY_MODES`].
pub const MAX_BODY_MODES: usize = 256;
pub const MAX_PLATE_BODY_MODES: usize = 16;
/// Modal coordinates of one complete instrument (strings, twist, body, cavity, gas memory). The shared exact
/// propagation holds at most 42 of them as one dense block; a larger network needs the explicit
/// `prepared-linear` record, whose own library ceiling is `MAX_PREPARED_MODAL_MODES` (1024). The network
/// setup work `modes * (links + 1)^2` still has to fit [`MAX_WORK`].
pub const MAX_TOTAL_MODES: usize = 1024;
const _: () = assert!(MAX_TOTAL_MODES <= fs_couple::render::schedule::force::coupled::MAX_PREPARED_MODAL_MODES);
pub const MAX_EVENTS: usize = 16_384;
pub const MAX_CONTROLS: usize = 262_144;
/// The complete network's setup work `modes * (links + 1)^2` at the ceilings: every mode of the total ceiling with one
/// link per body mode at the body ceiling. It also bounds plate and relaxation projection terms; each performance's
/// declared budgets bound the real work, and this ceiling only refuses impossible declarations.
pub const MAX_WORK: u64 = (MAX_TOTAL_MODES * (MAX_BODY_MODES + 1) * (MAX_BODY_MODES + 1)) as u64;

#[derive(Clone, Debug)]
pub struct Wire {
    pub pack: String, pub length: f64, pub radius: f64, pub prestress: StringPrestress,
    pub section: WireSection,
    pub bending_loss: BendingLoss,
    pub torsion_memory: Option<Vec<torsion_memory::LayerProperties>>,
    /// With fixed extension, L0 is stress-free at this explicit reference T.
    /// Bonded layers share that reference but retain their own expansion curves.
    /// None preserves the original same-temperature mechanical prescription.
    pub thermal_reference_kelvin: Option<f64>,
    pub bending_modes: usize, pub twist_modes: usize, pub bridge: SharedPlateBridge,
    pub polarizations: Option<TransverseDirections>,
    /// The `tension-loss` row this string selects. None is no tension loss term.
    pub tension_loss: Option<TensionLossRow>,
}
/// One cited construction's tension loss factor eta_F in `tension-loss.tsv`.
#[derive(Clone, Debug, PartialEq)]
pub struct TensionLossRow { pub id: String, pub eta_f: f64, pub source: String }
/// The evidence table that `tension-loss ID` selects from, compiled into the example.
pub const TENSION_LOSS_TABLE: &str = include_str!("tension-loss.tsv");
/// Every row of the table in file order. A malformed, duplicate or nonphysical row refuses.
pub fn tension_loss_rows() -> Result<Vec<TensionLossRow>> {
    let mut lines = TENSION_LOSS_TABLE.lines().filter(|line| !line.starts_with('#'));
    if lines.next() != Some("id\teta_F\tsource") {
        return Err("tension-loss.tsv needs the header id, eta_F, source".into());
    }
    let mut rows: Vec<TensionLossRow> = Vec::new();
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        let &[id, eta_f, source] = fields.as_slice() else {
            return Err("tension-loss.tsv rows need exactly id, eta_F and source".into());
        };
        let eta_f: f64 = eta_f.parse()?;
        if id.is_empty() || id.contains(char::is_whitespace) || source.is_empty() || !eta_f.is_finite()
            || eta_f < 0.0 || rows.iter().any(|row| row.id == id) {
            return Err("tension-loss.tsv rows need a unique id, a finite nonnegative eta_F and a source".into());
        }
        rows.push(TensionLossRow { id: id.into(), eta_f, source: source.into() });
    }
    Ok(rows)
}
/// A source-property selection; the original specimen owns every loss law.
#[derive(Clone, Debug)]
pub enum BendingLoss {
    KelvinVoigt,
    Prony { properties: Vec<BendingRelaxationProperties> },
}
impl BendingLoss {
    fn validate(&self) -> Result<()> {
        let Self::Prony { properties } = self else { return Ok(()); };
        if !(1..=8).contains(&properties.len()) {
            return Err("string Prony bending needs 1..=8 explicit modulus/time pairs".into());
        }
        let mut moduli = std::collections::BTreeSet::new();
        let mut times = std::collections::BTreeSet::new();
        for (i, pair) in properties.iter().enumerate() {
            let modulus = pair.modulus.as_str();
            let time = pair.relaxation_time.as_str();
            if modulus.is_empty() || time.is_empty() || modulus.len() > 256 || time.len() > 256
                || modulus.split_ascii_whitespace().count() != 1
                || time.split_ascii_whitespace().count() != 1
                || modulus == time || properties[..i].contains(pair)
                || times.contains(modulus) || moduli.contains(time) {
                return Err("string memory requires distinct bounded property pairs with consistent physical quantities".into());
            }
            moduli.insert(modulus); times.insert(time);
        }
        Ok(())
    }
}
/// Section construction is explicit: bonding and a sliding close-coiled helix
/// have different mass, extension and bending laws.
#[derive(Clone, Debug)]
pub enum WireSection {
    Homogeneous,
    /// Whole-assembly measurements with an explicit elastic-only assumption.
    /// Coordinates belong to the measurement, independently of body/gas state.
    Measured { axes: Vec<(String, f64)>, properties: MeasuredStringProperties },
    Bonded { core_radius_m: f64, sleeves: Vec<WireLayer> },
    CloseWound { core_radius_m: f64, pack: String, wire_radius_m: f64, pitch_m: f64 },
}
#[derive(Clone, Debug)]
pub struct WireLayer { pub pack: String, pub outer_radius_m: f64 }
const MAX_LAYERS: usize = 16;
const MAX_ROTATIONAL_BRIDGES: usize = 16;
const MAX_TRANSLATIONAL_BRIDGES: usize = 16;
/// An additional physical displacement axis at a node of the retained shell.
/// String directions use its transverse frame; body directions use global XYZ.
#[derive(Clone, Debug, PartialEq)]
pub struct TranslationalBridge {
    pub string: usize, pub station_fraction: f64, pub direction: TransverseDirection,
    pub shell_component: usize, pub node: usize, pub body_direction: [f64; 3],
    pub stiffness_n_m: f64, pub damping_n_s_m: f64, pub rest_extension_m: f64,
}
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

/// One geometry-derived exterior load on the original material body basis.
/// The source description chooses the physical boundary problem, never a name
/// for an instrument or a material-dependent output filter.
#[derive(Clone, Copy, Debug)]
pub enum BodyRadiation {
    Rayleigh(PlateRayleighOptions),
    ShellBem(ShellBemOptions),
    /// The same Galerkin load on the body profile's declared acoustic skin.
    Skin(SkinRadiation),
}
/// The two-way load on the declared skin: the body modes inside the band, then
/// an exterior-owned neck through the skin's aperture caps.
#[derive(Clone, Copy, Debug)]
pub struct SkinRadiation {
    pub options: ShellBemOptions,
    pub angular_refinement: u32,
    pub maximum_anchor_distance_m: f64,
}
impl BodyRadiation {
    pub fn validate(self, rate: u32) -> Result<()> {
        match self {
            Self::Rayleigh(options) => options.validate(rate)?,
            Self::ShellBem(options) => options.validate(rate)?,
            Self::Skin(load) => load.options.validate(rate)?,
        }
        Ok(())
    }
    pub fn maximum_poles(self) -> usize {
        match self {
            Self::Rayleigh(options) => options.maximum_poles,
            Self::ShellBem(options) => options.maximum_poles,
            Self::Skin(load) => load.options.maximum_poles,
        }
    }
}
/// Explicit local caloric reservoirs of the original cavity atmosphere.
#[derive(Clone, Debug)]
pub struct CavityMolecular {
    pub relaxation: RelaxationState,
    pub options: ShellCavityMolecularOptions,
}

#[derive(Clone, Debug)]
pub struct Spec {
    pub bytes: Vec<u8>, pub rate: u32, pub samples: u64, pub full_scale: f64, pub block: usize,
    pub temperature: f64, pub omega: f64, pub gas: GasState, pub body: Body, pub strings: Vec<Wire>,
    /// Additional held source coordinates, keyed by the original material or
    /// bow/bridge/axial interface pack path. Each physical owner retains its
    /// temperature, frequency, pressure and actual-slip conventions.
    pub material_states: BTreeMap<String, Vec<(String, f64)>>,
    /// Retain original thermochemical composition and pack authority with the coefficients.
    pub gas_record: GasRecord,
    /// Optional independent interior atmosphere. Omission inherits gas_record;
    /// an explicit record must be consumed by a lumped or spatial body cavity.
    pub cavity_gas_record: Option<GasRecord>,
    /// Interior molecular data are independent of exterior observer absorption.
    pub cavity_molecular: Option<CavityMolecular>,
    pub bridge_interfaces: bridge_interfaces::Inputs,
    pub rotational_bridges: Vec<RotationalBridge>,
    pub translational_bridges: Vec<TranslationalBridge>,
    pub cylinder_gas_loading: Option<CylinderFlowResolution>,
    pub prepared_linear: Option<PreparedModalBudget>,
    /// `residual-flexibility exact-static-green`: the omitted bending modes' static compliance in series with
    /// every bridge and fixed stop (pexfp.299). Absent means the retained modes alone.
    pub residual_flexibility: bool,
    pub geometric_stretching: Option<(f64, NonlinearModalConfig)>,
    pub longitudinal_waves: Option<(LongitudinalStringResolution, NonlinearModalConfig)>,
    pub longitudinal_loss: LongitudinalMaterialLoss,
    pub longitudinal_memory: Option<longitudinal_memory::Input>,
    pub axial_bearings: axial_bearings::Inputs,
    pub bending_memory: Option<bending_memory::Input>,
    pub body_radiation: Option<BodyRadiation>,
    pub shell_bem_angular_refinement: u32,
    /// The same loaded plate surface with material gas on every receiver path.
    pub plate_observer: Option<PlateObserverOptions>,
    pub shell_observer: Option<ShellObserverOptions>,
    pub shell_receiver_quadrature: fs_bem::near_field::Options,
    /// One-way receiver through the body profile's closed acoustic skin; nothing is fed back.
    pub one_way_observer: Option<OneWayShellObserverOptions>,
    pub cavity: Option<MaterialBodyCavity>, pub friction: bow_interfaces::Inputs,
    /// The explicit contact-thermal profile requires original film density.
    pub contact_thermal: bool,
    /// An explicit profile routes dry heat through original material-coordinate strips.
    pub dry_contact_thermal: bool,
    /// Every film has a v3 route to physical fibers; its mass/history stays on the string.
    pub retained_string_films: bool,
    pub bow_ribbon: Option<(f64, usize)>, pub forces: Option<MaterialForceProgram>,
    pub bow_fibers: Option<fibers::Input>,
    /// Prescribed position along the string, separate from tangential speed/load.
    pub bow_station_motion: Option<bow_motion::Input>,
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
        if rate == 0 || rate > 384_000 || samples == 0 || samples > 28_800_000
            || samples > u64::from(rate)*600 || block == 0 {
            return Err("audio requires a bounded positive clock, horizon and callback".into());
        }
        let mut row = r.row("state")?;
        let temperature = row.positive()?; let omega = row.scalar()?; row.end()?;
        if omega < 0.0 { return Err("source angular frequency must be nonnegative".into()); }
        let mut material_states = BTreeMap::new();
        if r.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next()) == Some("material-states") {
            let mut row = r.row("material-states")?;
            let count = row.count(64)?; row.end()?;
            for _ in 0..count {
                let mut row = r.row("material-state")?;
                let pack = row.path()?;
                let count = row.count(8)?;
                if count == 0 { return Err("a material-state record needs at least one coordinate".into()); }
                let mut axes = Vec::with_capacity(count);
                for _ in 0..count { axes.push((row.label()?, row.scalar()?)); }
                row.end()?;
                if material_states.insert(pack, axes).is_some() {
                    return Err("material-state records must name distinct pack paths".into());
                }
            }
        }
        let mut row = r.row("gas")?;
        let gas_record = read_gas_record(&mut row.words, None, sources)?; row.end()?;
        let gas = *gas_record.state();
        let cavity_gas_record = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("cavity_gas") {
            let mut row = r.row("cavity_gas")?;
            let record = read_gas_record(&mut row.words, None, sources)?;
            row.end()?;
            fs_couple::cavity_modal::validate_cavity_media(record.state(), &gas)?;
            Some(record)
        } else { None };
        let cavity_gas = *cavity_gas_record.as_ref().unwrap_or(&gas_record).state();
        let cavity_molecular = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("cavity_molecular") {
            let mut row = r.row("cavity_molecular")?;
            let options = ShellCavityMolecularOptions {
                maximum_caloric_coordinates: row.count(128)?,
                maximum_memory_states: row.count(4096)?,
            };
            let law = fs_material::relaxation::input::read_relaxation(&mut row.words)?;
            row.end()?;
            let medium = fs_material::relaxation::input::medium_from_record(
                cavity_gas_record.as_ref().unwrap_or(&gas_record))?;
            Some(CavityMolecular { relaxation: law.evaluate(&medium)?, options })
        } else { None };
        let mut row = r.row("body")?;
        let body = Body { pack: row.path()?, length: row.positive()?, width: row.positive()?,
            thickness: row.positive()?, modes: row.count(MAX_BODY_MODES)?, listener: row.positive()? };
        row.end()?;
        let mut row = r.row("strings")?; let count = row.count(8)?; row.end()?;
        if count == 0 || body.modes == 0 { return Err("nonempty material strings and body modes required".into()); }
        let mut strings = Vec::with_capacity(count);
        let mut bridge_sources = Vec::new();
        for string in 0..count {
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
            let bending_modes = row.count(MAX_MODES)?;
            let twist_modes = row.count(MAX_MODES)?;
            let station_fraction = row.scalar()?;
            let sourced = row.words.clone().next() == Some("interface");
            let (stiffness_n_m, damping_n_s_m) = if sourced {
                row.word()?;
                bridge_sources.push(bridge_interfaces::Declaration {
                    string, pack: row.path()?, area_m2: row.positive()?,
                    stiffness_property: row.word()?.to_string(), viscosity_property: row.word()?.to_string(),
                });
                // Unresolved source inputs cannot reach a physics constructor;
                // the original interface owner supplies both coefficients first.
                (0.0, 0.0)
            } else { (row.scalar()?, row.scalar()?) };
            let mut wire = Wire { pack, length, radius, prestress, section: WireSection::Homogeneous,
                bending_loss: BendingLoss::KelvinVoigt, torsion_memory: None,
                thermal_reference_kelvin: None, tension_loss: None,
                polarizations: None, bending_modes, twist_modes,
                bridge: SharedPlateBridge { station_fraction, stiffness_n_m, damping_n_s_m, rest_extension_m: 0.0 } };
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
            if r.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next()) == Some("bending-loss") {
                if matches!(&wire.section, WireSection::Measured { .. }) {
                    return Err("measured whole strings do not supply an intrinsic bending-loss law".into());
                }
                let mut row = r.row("bending-loss")?;
                wire.bending_loss = match row.word()? {
                    "kelvin-voigt" => BendingLoss::KelvinVoigt,
                    "prony" => {
                        let count = row.count(8)?;
                        let mut properties = Vec::with_capacity(count);
                        for _ in 0..count {
                            properties.push(BendingRelaxationProperties {
                                modulus: row.word()?.into(), relaxation_time: row.word()?.into(),
                            });
                        }
                        BendingLoss::Prony { properties }
                    },
                    _ => return Err("bending-loss must be kelvin-voigt or prony COUNT MODULUS TIME ...".into()),
                };
                row.end()?;
                wire.bending_loss.validate()?;
            }
            if r.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next()) == Some("tension-loss") {
                let mut row = r.row("tension-loss")?;
                let id = row.word()?;
                row.end()?;
                wire.tension_loss = Some(tension_loss_rows()?.into_iter().find(|row| row.id == id)
                    .ok_or("tension-loss names no row of tension-loss.tsv")?);
            }
            if r.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next()) == Some("torsion-memory") {
                let mut row = r.row("torsion-memory")?;
                let count = row.count(MAX_LAYERS)?; row.end()?;
                let mut layers = Vec::with_capacity(count);
                for _ in 0..count {
                    let mut row = r.row("torsion-layer")?;
                    let equilibrium_modulus = row.word()?.to_string();
                    let count = row.count(8)?;
                    let mut branches = Vec::with_capacity(count);
                    for _ in 0..count {
                        branches.push((row.word()?.to_string(), row.word()?.to_string()));
                    }
                    row.end()?;
                    layers.push(torsion_memory::LayerProperties { equilibrium_modulus, branches });
                }
                wire.torsion_memory = Some(layers);
                wire.torsion_memory_states()?;
            }
            if wire.bending_modes == 0 || !(0.0 < wire.bridge.station_fraction && wire.bridge.station_fraction < 1.0)
                || wire.bridge.stiffness_n_m < 0.0 || wire.bridge.damping_n_s_m < 0.0
                || (!sourced && !(wire.bridge.stiffness_n_m > 0.0 || wire.bridge.damping_n_s_m > 0.0)) {
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
                    || !matches!(wire.prestress, StringPrestress::FixedExtension { .. }) {
                    return Err("thermal support requires one declaration per fixed-extension string".into());
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
        let mut translational_bridges = Vec::new();
        if r.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next()) == Some("translational-bridges") {
            let mut row = r.row("translational-bridges")?;
            let joints = row.count(MAX_TRANSLATIONAL_BRIDGES)?; row.end()?;
            for _ in 0..joints {
                let mut row = r.row("translational-bridge")?;
                let joint = TranslationalBridge {
                    string: row.count(count - 1)?, station_fraction: row.scalar()?,
                    direction: TransverseDirection::try_new([row.scalar()?, row.scalar()?])?,
                    shell_component: row.count(63)?, node: row.count(3071)?,
                    body_direction: [row.scalar()?, row.scalar()?, row.scalar()?],
                    stiffness_n_m: row.scalar()?, damping_n_s_m: row.scalar()?,
                    rest_extension_m: row.scalar()?,
                };
                row.end()?;
                let norm2 = joint.body_direction.iter().map(|x| x*x).sum::<f64>();
                if !(0.0 < joint.station_fraction && joint.station_fraction < 1.0)
                    || !norm2.is_finite() || (norm2 - 1.0).abs() > 1e-12
                    || (strings[joint.string].polarizations.is_none() && joint.direction.components()[1] != 0.0)
                    || joint.stiffness_n_m < 0.0 || joint.damping_n_s_m < 0.0
                    || (joint.stiffness_n_m == 0.0 && joint.damping_n_s_m == 0.0) {
                    return Err("translational bridge needs an interior station, retained transverse direction, unit global body axis and a nonzero passive law".into());
                }
                translational_bridges.push(joint);
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
        let residual_flexibility = r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("residual-flexibility");
        if residual_flexibility {
            let mut row = r.row("residual-flexibility")?;
            if row.word()? != "exact-static-green" {
                return Err("residual-flexibility names its one source: exact-static-green".into());
            }
            row.end()?;
        }
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
        let longitudinal_waves = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("longitudinal-waves") {
            let mut row = r.row("longitudinal-waves")?;
            let maximum_slope = row.positive()?;
            let maximum_strain = row.positive()?;
            let config = NonlinearModalConfig {
                max_iterations: row.count(128)?, maximum_generalized_force: row.positive()?,
                force_absolute_tolerance: row.positive()?, force_relative_tolerance: row.positive()?,
            };
            let count = row.count(32)?;
            let mut axial_modes = Vec::with_capacity(count);
            for _ in 0..count { axial_modes.push(row.count(1024)?); }
            row.end()?;
            let resolution = LongitudinalStringResolution { axial_modes, maximum_slope, maximum_strain };
            resolution.validate()?;
            config.validate()?;
            if geometric_stretching.is_some() {
                return Err("longitudinal waves replace geometric stretching; select one strain law".into());
            }
            Some((resolution, config))
        } else { None };
        let longitudinal_loss = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("longitudinal-loss") {
            if longitudinal_waves.is_none() {
                return Err("longitudinal-loss requires explicit longitudinal-waves".into());
            }
            let mut row = r.row("longitudinal-loss")?;
            let loss = match row.word()? {
                "elastic" => LongitudinalMaterialLoss::Elastic,
                "kelvin-voigt" => LongitudinalMaterialLoss::KelvinVoigt {
                    viscosity_property: fs_material::section::axial::YOUNG_VISCOSITY_PROPERTY.into(),
                },
                "finite-strain-kelvin-voigt" => LongitudinalMaterialLoss::FiniteStrainKelvinVoigt {
                    viscosity_property: fs_material::section::axial::YOUNG_VISCOSITY_PROPERTY.into(),
                },
                _ => return Err("longitudinal-loss must be elastic, kelvin-voigt or finite-strain-kelvin-voigt".into()),
            };
            row.end()?;
            loss
        } else { LongitudinalMaterialLoss::Elastic };
        let longitudinal_memory = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("longitudinal-memory") {
            let mut row = r.row("longitudinal-memory")?;
            let finite_strain=row.words.clone().next()==Some("finite-strain");
            if finite_strain {row.word()?;}
            let max_states = row.count(256)?;
            let equilibrium_modulus = row.word()?.to_string();
            let count = row.count(8)?;
            let mut branches = Vec::with_capacity(count);
            for _ in 0..count {
                branches.push((row.word()?.to_string(), row.word()?.to_string()));
            }
            row.end()?;
            let memory = longitudinal_memory::Input { max_states, finite_strain,
                properties: fs_material::section::relaxation::LayerRelaxationProperties {
                    equilibrium_modulus, branches,
                },
            };
            memory.validate()?;
            Some(memory)
        } else { None };
        let mut axial_declarations = Vec::new();
        if r.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next()) == Some("axial-bearings") {
            use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::longitudinal::{
                AxialBearingGeometry, bearings::MAX_AXIAL_BEARINGS,
            };
            if longitudinal_waves.is_none() {
                return Err("axial bearings require explicit longitudinal-waves".into());
            }
            let mut row = r.row("axial-bearings")?;
            let bearings = row.count(MAX_AXIAL_BEARINGS)?; row.end()?;
            for _ in 0..bearings {
                let mut row = r.row("axial-bearing")?;
                let geometry = AxialBearingGeometry {
                    string: row.count(count-1)?, station_fraction: row.scalar()?,
                    shell_component: row.count(63)?, node: row.count(3071)?,
                    string_axis: [row.scalar()?, row.scalar()?, row.scalar()?],
                    contact_arm_m: [row.scalar()?, row.scalar()?, row.scalar()?],
                    rest_extension_m: row.scalar()?,
                };
                geometry.validate()?;
                if row.word()? != "interface" {
                    return Err("axial bearing requires an explicit sourced interface".into());
                }
                let source = bridge_interfaces::Declaration {
                    string: geometry.string, pack: row.path()?, area_m2: row.positive()?,
                    stiffness_property: row.word()?.into(), viscosity_property: row.word()?.into(),
                };
                row.end()?;
                axial_declarations.push(axial_bearings::Declaration { geometry, source });
            }
        }
        let axial_bearings = axial_bearings::Inputs::new(axial_declarations)?;
        let bending_memory = if r.lines.clone().next()
            .and_then(|line|line.split_ascii_whitespace().next()) == Some("bending-memory") {
            let mut row=r.row("bending-memory")?;
            let max_states=row.count(256)?; let selected=row.count(count)?; row.end()?;
            let mut selections=Vec::with_capacity(selected);
            for _ in 0..selected {
                let mut row=r.raw_row()?;
                match row.word()? {
                    "string-memory" => {
                        let string=row.count(count-1)?; let terms=row.count(8)?;
                        let mut properties=Vec::with_capacity(terms);
                        for _ in 0..terms { properties.push((row.word()?.to_string(),row.word()?.to_string())); }
                        row.end()?;
                        selections.push(bending_memory::Input::uniform_selection(string,properties,&strings)?);
                    },
                    "layered-string-memory" => {
                        let string=row.count(count-1)?; let count=row.count(MAX_LAYERS)?;
                        row.end()?;
                        let mut layers=Vec::with_capacity(count);
                        for _ in 0..count {
                            let mut row=r.row("layer-memory")?;
                            let equilibrium_modulus=row.word()?.to_string(); let terms=row.count(8)?;
                            let mut branches=Vec::with_capacity(terms);
                            for _ in 0..terms { branches.push((row.word()?.to_string(),row.word()?.to_string())); }
                            row.end()?;
                            layers.push(bending_memory::LayerProperties {equilibrium_modulus,branches});
                        }
                        selections.push(bending_memory::Selection {string,layers});
                    },
                    _ => return Err("expected string-memory or layered-string-memory".into()),
                }
            }
            Some(bending_memory::Input::new(max_states,selections,&strings)?)
        } else { None };
        let body_radiation = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("body-radiation") {
            let mut row = r.row("body-radiation")?;
            let options = match row.word()? {
                "rayleigh" => BodyRadiation::Rayleigh(PlateRayleighOptions {
                    sides: row.parse()?, minimum_hz: row.positive()?, maximum_hz: row.positive()?,
                    frequency_samples: row.count(257)?, maximum_surface_cells: row.count(4096)?,
                    maximum_sweep_terms: row.parse()?, maximum_poles: row.count(256)?,
                }),
                "bem" => BodyRadiation::ShellBem(ShellBemOptions {
                    minimum_hz: row.positive()?, maximum_hz: row.positive()?,
                    frequency_samples: row.count(257)?, maximum_surface_panels: row.count(2048)?,
                    maximum_sweep_terms: row.parse()?, maximum_poles: row.count(256)?,
                }),
                "skin" => BodyRadiation::Skin(SkinRadiation {
                    options: ShellBemOptions {
                        minimum_hz: row.positive()?, maximum_hz: row.positive()?,
                        frequency_samples: row.count(257)?,
                        maximum_surface_panels: row.count(fs_bem::reciprocal_radiation::MAX_RECIPROCAL_PANELS)?,
                        maximum_sweep_terms: row.parse()?,
                        maximum_poles: row.count(fs_couple::render::plate::impact::radiation::MAX_POLES)?,
                    },
                    angular_refinement: row.count(fs_bem::reciprocal_radiation::MAX_ANGULAR_REFINEMENT as usize)? as u32,
                    maximum_anchor_distance_m: row.positive()?,
                }),
                _ => return Err("body-radiation requires rayleigh, bem or skin geometry-derived loading".into()),
            };
            row.end()?;
            options.validate(rate)?;
            // Construction selects ordinary/prepared playback or harmonic-only
            // analysis after parsing. Keep the physical source declaration;
            // that selected owner checks backend compatibility before pack I/O.
            Some(options)
        } else { None };
        let shell_bem_angular_refinement = if r.lines.clone().next()
            .and_then(|line|line.split_ascii_whitespace().next())==Some("shell_angular_refinement") {
            if !matches!(body_radiation,Some(BodyRadiation::ShellBem(_))) {
                return Err("shell_angular_refinement requires body-radiation bem".into());
            }
            let mut row=r.row("shell_angular_refinement")?;
            let level=row.count(2)? as u32;row.end()?;level
        } else {0};
        let plate_observer = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("plate_observer") {
            let mut row = r.row("plate_observer")?;
            let position_m = [row.scalar()?, row.scalar()?, row.scalar()?];
            let maximum_frequency_hz = row.positive()?;
            let maximum_surface_cells = row.count(4096)?;
            let maximum_replay_terms = row.count(fs_vfit::observation::MAX_OBSERVATION_ENTRIES)?;
            row.end()?;
            let mut row = r.row("observer_gas")?;
            let propagation = fs_couple::pcm_wav::baffled::gas::GasPropagation::read_record_for_gas(
                &mut row.words, &gas_record)?;
            let options = PlateObserverOptions { position_m, maximum_frequency_hz,
                maximum_surface_cells, maximum_replay_terms, propagation };
            options.validate(rate)?;
            Some(options)
        } else { None };
        let shell_observer = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("shell_observer") {
            let mut row = r.row("shell_observer")?;
            let receiver = fs_bem::observation_transfers::ReceiverSpec {
                id: row.word()?.to_string(),
                position_m: [row.scalar()?, row.scalar()?, row.scalar()?],
            };
            let minimum_clearance_m = row.positive()?;
            let pass_band_upper_hz = row.positive()?;
            let taper_stop_hz = row.positive()?;
            let filter_latency_samples =
                row.count(fs_bem::observation_transfers::MAX_OBSERVATION_FFT / 4)?;
            let count = row.count(fs_bem::observation_transfers::MAX_OBSERVATION_RUNGS)?;
            let mut fft_lengths = Vec::with_capacity(count);
            for _ in 0..count {
                fft_lengths.push(row.count(fs_bem::observation_transfers::MAX_OBSERVATION_FFT)?);
            }
            let tail_energy_tolerance = row.positive()?;
            let dtft_relative_tolerance = row.positive()?;
            let maximum_frequency_entries =
                row.count(fs_vfit::observation::MAX_OBSERVATION_ENTRIES)?;
            let maximum_sweep_terms = row.parse()?;
            let maximum_replay_terms = row.parse()?;
            row.end()?;
            Some(ShellObserverOptions {
                receiver, minimum_clearance_m, pass_band_upper_hz, taper_stop_hz,
                filter_latency_samples, fft_lengths, tail_energy_tolerance,
                dtft_relative_tolerance, maximum_frequency_entries,
                maximum_sweep_terms, maximum_replay_terms,
            })
        } else { None };
        let shell_receiver_quadrature = if r.lines.clone().next()
            .and_then(|line|line.split_ascii_whitespace().next())==Some("shell_quadrature") {
            if shell_observer.is_none() { return Err("shell_quadrature requires shell_observer".into()); }
            let mut row=r.row("shell_quadrature")?;
            let options=fs_bem::near_field::Options {
                relative_tolerance:row.positive()?,maximum_depth:row.count(16)?,
                maximum_kernel_evaluations:row.count(20_000_000)?,
            };
            row.end()?;options.validate()?;options
        } else {fs_bem::near_field::Options::default()};
        let one_way_observer = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("one_way_observer") {
            let mut row = r.row("one_way_observer")?;
            let receiver = fs_bem::observation_transfers::ReceiverSpec {
                id: row.word()?.to_string(),
                position_m: [row.scalar()?, row.scalar()?, row.scalar()?],
            };
            let minimum_clearance_m = row.positive()?;
            let pass_band_upper_hz = row.positive()?;
            let taper_stop_hz = row.positive()?;
            let filter_latency_samples =
                row.count(fs_bem::observation_transfers::MAX_OBSERVATION_FFT / 4)?;
            let count = row.count(fs_bem::observation_transfers::MAX_OBSERVATION_RUNGS)?;
            let mut fft_lengths = Vec::with_capacity(count);
            for _ in 0..count {
                fft_lengths.push(row.count(fs_bem::observation_transfers::MAX_OBSERVATION_FFT)?);
            }
            let tail_energy_tolerance = row.positive()?;
            let dtft_relative_tolerance = row.positive()?;
            let maximum_frequency_entries =
                row.count(fs_vfit::observation::MAX_OBSERVATION_ENTRIES)?;
            let maximum_replay_terms = row.parse()?;
            let panel_rule = match row.word()? {
                "triangle-near" => fs_bem::helmholtz::PanelRule::TriangleNear { max_pairs: row.parse()? },
                "triangle" => fs_bem::helmholtz::PanelRule::Triangle { max_pairs: row.parse()? },
                _ => return Err("one_way_observer panel rule must be triangle-near MAX_PAIRS or triangle MAX_PAIRS".into()),
            };
            let maximum_skin_panels = row.count(fs_bem::helmholtz::MAX_DENSE_PANELS)?;
            let maximum_anchor_distance_m = row.positive()?;
            let maximum_boundary_terms = row.parse()?;
            let workers = row.count(256)?;
            let rule = OneWayRule { maximum_loss_fraction: row.scalar()?, maximum_shift_fraction: row.scalar()? };
            row.end()?;
            Some(OneWayShellObserverOptions {
                receiver, minimum_clearance_m, pass_band_upper_hz, taper_stop_hz,
                filter_latency_samples, fft_lengths, tail_energy_tolerance,
                dtft_relative_tolerance, maximum_frequency_entries, maximum_replay_terms,
                panel_rule, maximum_skin_panels, maximum_anchor_distance_m,
                maximum_boundary_terms, workers, rule,
            })
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
                Some(MaterialBodyCavity { cavity: HelmholtzModalCavity::new(description, cavity_gas,
                    CavityNeckLoss::Resistance(resistance))?, modal_budget: modal_budget(), listener_m })
            },
            kind @ ("physical-neck" | "exterior-neck") => {
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
                // Interior gas supplies inertia, compression and Bessel neck
                // loss. A compact outer-end law still requires that same gas
                // outside; a skin load can own a distinct exterior atmosphere.
                // The original owner counts each end correction exactly once.
                // An exterior-owned outer end gives its mass and radiation to the
                // skin's computed load on the aperture caps, and has none without it.
                let loss = if kind == "exterior-neck" {
                    if !matches!(body_radiation, Some(BodyRadiation::Skin(_))) {
                        return Err("exterior-neck has no outer mass without the body-radiation skin load on its apertures".into());
                    }
                    CavityNeckLoss::ExteriorMouth
                } else {
                    CavityNeckLoss::RadiationAtResonance
                };
                Some(MaterialBodyCavity { cavity: HelmholtzModalCavity::new(description, cavity_gas, loss)?,
                    modal_budget: modal_budget(), listener_m })
            },
            _ => return Err("cavity must be none, neck V R L RESISTANCE LISTENER_M, or physical-neck|exterior-neck V R L INNER_END OUTER_END LISTENER_M".into()),
        };
        row.end()?;
        if let Some(cavity) = &cavity { cavity.cavity.validate_exterior_gas(&gas)?; }
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
            "hyperbolic" => {
                let law=FrictionLaw::Hyperbolic { static_mu:row.scalar()?,kinetic_mu:row.scalar()?,
                    characteristic_speed:row.positive()? };
                row.end()?; bow_interfaces::Inputs::research(law)?
            },
            "layers" => {
                let count=row.count(16)?; row.end()?;
                let mut declarations=Vec::with_capacity(count);
                for _ in 0..count {
                    let mut row=r.row("shear-layer")?;
                    let mut declaration=shear_layers::Declaration {
                        track:row.word()?.to_string(), pack:row.path()?,
                        area_m2:row.positive()?, thickness_m:row.positive()?,
                        maximum_elastic_shear_strain:row.positive()?,
                        relaxation: Vec::new(),
                    };
                    row.end()?;
                    if r.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next())
                        == Some("shear-memory") {
                        let mut row = r.row("shear-memory")?;
                        let arms = row.count(fs_material::state_point::MAX_SHEAR_RELAXATION_ARMS)?;
                        if arms == 0 { return Err("shear-memory requires at least one source arm".into()); }
                        for _ in 0..arms {
                            declaration.relaxation.push((row.word()?.to_string(), row.word()?.to_string()));
                        }
                        row.end()?;
                    }
                    declarations.push(declaration);
                }
                bow_interfaces::Inputs::layers(declarations)?
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
                        "hyperbolic"=>bow_interfaces::Model::Hyperbolic,
                        _=>return Err("bow interface model must be coulomb, stribeck or hyperbolic".into()),
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
        let bow_fibers = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("bow_fibers") {
            let mut row = r.row("bow_fibers")?;
            let pack = row.path()?;
            let radius_m = row.positive()?;
            let length_m = row.positive()?;
            let fiber_count = row.count(1 << 24)?;
            let (prestress, thermal_reference_kelvin) = match row.words.clone().next() {
                Some("fixed-extension") | Some("thermal-extension") => {
                    let thermal = row.word()? == "thermal-extension";
                    let support = StringPrestress::FixedExtension {
                        stress_free_length_m: row.positive()?,
                        linear_strain_limit: row.positive()?,
                    };
                    let reference = if thermal { Some(row.positive()?) } else { None };
                    (support, reference)
                }
                _ => (StringPrestress::FixedTension(row.positive()?), None),
            };
            let mut input = fibers::Input { pack, radius_m, length_m, fiber_count,
                prestress, thermal_reference_kelvin, station_fraction: row.scalar()?,
                axial_modes: row.count(MAX_MODES)?, transverse_modes: row.count(MAX_MODES)?,
                memory: None, strain: None, motion: None, support_program: None, normal_program: None,
            };
            row.end()?;
            if r.lines.clone().next()
                .and_then(|line| line.split_ascii_whitespace().next()) == Some("bow_fiber_memory") {
                let mut row = r.row("bow_fiber_memory")?;
                let max_states = row.count(256)?;
                let equilibrium_modulus = row.word()?.to_string();
                let count = row.count(fs_material::section::relaxation::MAX_SECTION_RELAXATION_BRANCHES)?;
                let mut branches = Vec::with_capacity(count);
                for _ in 0..count {
                    branches.push((row.word()?.to_string(), row.word()?.to_string()));
                }
                row.end()?;
                input.memory = Some(fibers::Memory {
                    properties: fs_material::section::relaxation::LayerRelaxationProperties {
                        equilibrium_modulus, branches,
                    },
                    max_states,
                });
            }
            if r.lines.clone().next()
                .and_then(|line| line.split_ascii_whitespace().next()) == Some("bow_fiber_strain") {
                let mut row = r.row("bow_fiber_strain")?;
                let limits = fibers::FiberStrainLimits {
                    maximum_slope: row.positive()?, maximum_strain: row.positive()?,
                };
                let config = NonlinearModalConfig {
                    max_iterations: row.count(128)?, maximum_generalized_force: row.positive()?,
                    force_absolute_tolerance: row.positive()?, force_relative_tolerance: row.positive()?,
                };
                row.end()?;
                limits.validate()?;
                config.validate()?;
                input.strain = Some((limits, config));
            }
            if r.lines.clone().next()
                .and_then(|line| line.split_ascii_whitespace().next()) == Some("bow_fiber_motion") {
                let mut row = r.row("bow_fiber_motion")?;
                let limits = fibers::FiberSupportMotionLimits {
                    maximum_speed_m_s: row.positive()?, maximum_acceleration_m_s2: row.positive()?,
                    maximum_step_m: row.positive()?, max_update_terms: row.count(1_048_576)?,
                };
                row.end()?;
                limits.validate()?;
                input.motion = Some(limits);
            }
            for (record, endpoint_record, normal) in [
                ("bow_fiber_program", "fiber_speed", false),
                ("bow_fiber_normal_program", "fiber_normal_speed", true),
            ] {
                if r.lines.clone().next()
                    .and_then(|line| line.split_ascii_whitespace().next()) != Some(record) {
                    continue;
                }
                let mut row = r.row(record)?;
                let initial_velocity_m_s = row.scalar()?;
                let maximum_speed_m_s = row.positive()?;
                let maximum_acceleration_m_s2 = row.positive()?;
                let maximum_abs_translation_m = row.positive()?;
                let max_update_terms = row.count(1_048_576)?;
                let count = row.count(MAX_EVENTS)?;
                row.end()?;
                let mut knots = Vec::with_capacity(count);
                for _ in 0..count {
                    let mut row = r.row(endpoint_record)?;
                    let sample = row.word()?.parse::<u64>()?;
                    let velocity_m_s = row.scalar()?;
                    row.end()?;
                    knots.push(fibers::FiberSupportKnot { sample, velocity_m_s });
                }
                let program = fibers::FiberSupportProgram {
                    initial_velocity_m_s, knots, maximum_speed_m_s,
                    maximum_acceleration_m_s2, maximum_abs_translation_m,
                    max_knots: MAX_EVENTS, max_update_terms,
                };
                if normal { input.normal_program = Some(program); }
                else { input.support_program = Some(program); }
            }
            Some(input)
        } else { None };
        let bow_station_motion = bow_motion::Input::parse(&mut r)?;
        let forces = forces::parse(&mut r, count, |string| strings[string].polarizations.is_some(),
            rate, samples, block)?;
        if r.lines.next().is_some() { return Err("unexpected header record".into()); }
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
                    return Err("gesture events change bow speed/load; use bow_station_motion for position along the string".into());
                }
            }
        }
        friction.validate_tracks(&schedule)?;
        let spec = Self { bytes: bytes.to_vec(), rate, samples, full_scale, block, temperature, omega,
            gas, gas_record, cavity_gas_record, cavity_molecular, body, strings, material_states, bridge_interfaces: bridge_interfaces::Inputs::new(bridge_sources)?,
            rotational_bridges, translational_bridges, cylinder_gas_loading, prepared_linear, residual_flexibility, geometric_stretching,
            longitudinal_waves, longitudinal_loss, longitudinal_memory, axial_bearings, bending_memory, body_radiation, shell_bem_angular_refinement,
            plate_observer, shell_observer, shell_receiver_quadrature, one_way_observer, cavity, friction, contact_thermal: false, dry_contact_thermal: false, retained_string_films: false, bow_ribbon, bow_fibers,
            bow_station_motion, forces, schedule };
        spec.validate_material_states()?;
        spec.bow_station_program(&CancelGate::new_clock_free())?;
        spec.validate_cavity_gas()?;
        spec.validate_shell_observer()?;
        if let Some(layers) = spec.friction.shear_layers() { layers.preflight(&spec)?; }
        spec.network_dimensions(false)?;
        Ok(spec)
    }

    /// Bound the extra fixed environment independently of source-pack I/O.
    pub fn validate_material_states(&self) -> Result<()> {
        if self.material_states.len() > 64 {
            return Err("material states exceed 64 source packs".into());
        }
        for (pack, axes) in &self.material_states {
            if pack.is_empty() || pack.len() > 512 || !std::path::Path::new(pack).components()
                .all(|part| matches!(part, std::path::Component::Normal(_)))
                || axes.is_empty() || axes.len() > 8
            {
                return Err("material states require relative pack paths and one to eight coordinates".into());
            }
            for (index, (axis, value)) in axes.iter().enumerate() {
                if axis.is_empty() || axis.len() > 256 || axis.chars().any(char::is_control)
                    || axis == "T" || axis == "omega" || !value.is_finite()
                    || axes[..index].iter().any(|(prior, _)| prior == axis)
                {
                    return Err("held material coordinates must be distinct finite values outside reserved T/omega".into());
                }
            }
            if self.strings.iter().any(|wire| wire.pack == *pack
                && matches!(&wire.section, WireSection::Measured { .. }))
            {
                return Err("measured whole strings retain their own explicit observation coordinates; material-state cannot replace them".into());
            }
        }
        Ok(())
    }

    /// The complete source record for interior pressure, inertia and heat loss.
    /// It may share data with the exterior by inheritance, never by replacing
    /// an explicit original regional declaration.
    pub fn cavity_gas_record(&self) -> &GasRecord {
        self.cavity_gas_record.as_ref().unwrap_or(&self.gas_record)
    }

    pub fn cavity_gas(&self) -> &GasState { self.cavity_gas_record().state() }

    pub fn validate_cavity_gas(&self) -> Result<()> {
        fs_couple::cavity_modal::validate_cavity_media(self.cavity_gas(), &self.gas)?;
        if let Some(molecular) = &self.cavity_molecular {
            let width = 1 + molecular.relaxation.processes().iter()
                .filter(|process| process.heat_capacity_j_kg_k > 0.0).count();
            if molecular.relaxation.gas() != self.cavity_gas()
                || !(1..=128).contains(&molecular.options.maximum_caloric_coordinates)
                || !(1..=4096).contains(&molecular.options.maximum_memory_states)
                || width > molecular.options.maximum_caloric_coordinates {
                return Err("cavity_molecular requires its original interior gas and complete positive caloric/memory allowances".into());
            }
            if self.cavity.is_some() || self.prepared_linear.is_some() {
                return Err("cavity_molecular requires an exact spatial thermal cavity, without a separate lumped cavity or prepared-linear owner".into());
            }
        }
        if let Some(cavity) = &self.cavity {
            if cavity.cavity.gas() != self.cavity_gas() {
                return Err("body cavity must retain its original interior gas source record".into());
            }
            cavity.cavity.validate_exterior_gas(&self.gas)?;
        }
        Ok(())
    }

    /// The optional listener must observe the same explicitly loaded shell.
    /// Repeat this cheap source admission before material I/O for programmatic
    /// edits to an already parsed specification.
    pub fn validate_shell_observer(&self) -> Result<()> {
        if self.shell_bem_angular_refinement>2
            || (self.shell_bem_angular_refinement!=0
                && !matches!(self.body_radiation,Some(BodyRadiation::ShellBem(_)))) {
            return Err("shell angular refinement requires the original BEM owner and level 0..2".into());
        }
        self.shell_receiver_quadrature.validate()?;
        if let Some(observer) = &self.plate_observer {
            observer.validate(self.rate)?;
            if !matches!(self.body_radiation, Some(BodyRadiation::Rayleigh(_)))
                || self.cavity.is_some() || self.shell_observer.is_some() || self.one_way_observer.is_some() {
                return Err("plate_observer requires its original body-radiation rayleigh load, no cavity mouth, and no other observer".into());
            }
        }
        if let Some(observer) = &self.shell_observer {
            observer.validate(self.rate)?;
            let Some(BodyRadiation::ShellBem(load)) = self.body_radiation else {
                return Err("shell_observer requires the original body-radiation bem load".into());
            };
            let finest = *observer.fft_lengths.last()
                .ok_or("shell_observer requires its explicit FFT ladder")?;
            let first_positive_hz = f64::from(self.rate) / finest as f64;
            if self.cavity.is_some() || observer.taper_stop_hz > load.maximum_hz
                || first_positive_hz < load.minimum_hz {
                return Err("shell_observer requires no cavity and its positive-frequency band covered by the original BEM load".into());
            }
        }
        if let Some(observer) = &self.one_way_observer {
            observer.validate(self.rate)?;
            // The skin's own two-way load is observed through the same skin.
            if matches!(self.body_radiation, Some(BodyRadiation::Rayleigh(_) | BodyRadiation::ShellBem(_)))
                || self.shell_observer.is_some() {
                return Err("one_way_observer takes only its own skin's computed load, neither another air load nor a shell_observer".into());
            }
        }
        Ok(())
    }

    pub fn bow_law(&self,port:usize)->Result<&FrictionLaw> {
        let track=self.schedule.tracks().get(port).ok_or("unknown bow contact")?;
        self.friction.law(&track.id)
    }

    pub fn bow_contact_count(&self)->usize {
        self.bow_ribbon.map_or(self.schedule.tracks().len(), |(_, cells)| cells)
    }

    /// Reserve every declared original-specimen arm before material I/O.
    /// Explicit zero strengths may need fewer states in the physical owner.
    pub fn string_memory_states(&self) -> Result<usize> {
        let mut states = 0usize;
        for wire in &self.strings {
            wire.bending_loss.validate()?;
            if let BendingLoss::Prony { properties } = &wire.bending_loss {
                let layers = match &wire.section {
                    WireSection::Bonded { sleeves, .. } => sleeves.len() + 1,
                    WireSection::Measured { .. } => return Err("measured whole strings do not supply layer relaxation laws".into()),
                    WireSection::Homogeneous | WireSection::CloseWound { .. } => 1,
                };
                if self.prepared_linear.is_some() {
                    return Err("string material memory requires the exact owner; prepared-linear has no internal-strain states".into());
                }
                let arms = layers.checked_mul(properties.len())
                    .ok_or("complete layer relaxation-arm count overflow")?;
                if arms > fs_material::section::relaxation::MAX_SECTION_RELAXATION_BRANCHES {
                    return Err("bonded string exceeds the original section relaxation-arm ceiling".into());
                }
                let planes = if wire.polarizations.is_some() { 2 } else { 1 };
                states = wire.bending_modes.checked_mul(planes)
                    .and_then(|n| n.checked_mul(arms))
                    .and_then(|n| states.checked_add(n))
                    .ok_or("complete string memory-state count overflow")?;
                if states > 256 {
                    return Err("complete declared string memory exceeds the original exact-state ceiling".into());
                }
            }
            if wire.torsion_memory.is_some() && self.prepared_linear.is_some() {
                return Err("torsional material memory requires the exact owner; prepared-linear has no internal-strain states".into());
            }
            states = states.checked_add(wire.torsion_memory_states()?)
                .ok_or("complete string memory-state count overflow")?;
            if states > 256 {
                return Err("complete declared string memory exceeds the original exact-state ceiling".into());
            }
        }
        Ok(states)
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
        self.network_dimensions_with_body_modes(tensor_body, additional_body_connections, 0)
    }

    /// Include every source-declared interior-air coordinate and compression
    /// column in the same complete mechanical allocation. Uniform compression
    /// contributes a connection but no fictitious inertia.
    pub fn network_dimensions_with_body_modes(&self, tensor_body: bool,
        additional_body_connections: usize, additional_body_modes: usize)
        -> Result<(usize, usize, usize)>
    {
        self.string_memory_states()?;
        if let Some(memory)=&self.bending_memory { memory.preflight(self)?; }
        if let Some(memory)=&self.longitudinal_memory { memory.preflight(self)?; }
        let mut modes = self.body.modes.checked_add(usize::from(self.cavity.is_some()))
            .and_then(|n| n.checked_add(additional_body_modes))
            .ok_or("complete mode count overflow")?;
        if let Some(fibers) = &self.bow_fibers {
            modes = modes.checked_add(fibers.mode_count(self)?).ok_or("complete fiber mode count overflow")?;
        }
        let mut fluid_loaded_string_modes = 0usize;
        for wire in &self.strings {
            let planes = if wire.polarizations.is_some() { 2 } else { 1 };
            let bending = wire.bending_modes.checked_mul(planes)
                .ok_or("transverse mode count overflow")?;
            // Every transverse and torsional coordinate has its own causal
            // cylinder-fluid memory. Count both before enforcing exact budgets.
            fluid_loaded_string_modes = fluid_loaded_string_modes.checked_add(bending)
                .and_then(|n| n.checked_add(wire.twist_modes))
                .ok_or("fluid-loaded string mode count overflow")?;
            modes = modes.checked_add(bending).and_then(|n| n.checked_add(wire.twist_modes))
                .ok_or("complete string/body mode count overflow")?;
        }
        let cavity_links = self.cavity.as_ref().map_or(0, |c|
            1 + usize::from(c.cavity.neck_resistance_n_s_m() != 0.0));
        let mut links = self.strings.len().checked_add(self.rotational_bridges.len())
            .and_then(|n| n.checked_add(self.translational_bridges.len()))
            .and_then(|n| n.checked_add(self.axial_bearings.len()))
            .and_then(|n| n.checked_add(cavity_links))
            .and_then(|n| n.checked_add(if tensor_body { self.body.modes } else { 0 }))
            .and_then(|n| n.checked_add(additional_body_connections))
            .ok_or("complete connection count overflow")?;
        if let Some(resolution) = self.cylinder_gas_loading {
            resolution.validate()?;
            let memory = fluid_loaded_string_modes.checked_mul(resolution.intervals)
                .ok_or("fluid memory dimension overflow")?;
            modes = modes.checked_add(memory).ok_or("complete fluid mode count overflow")?;
            links = links.checked_add(memory).ok_or("complete fluid connection count overflow")?;
        }
        if self.axial_bearings.len() != 0 && self.longitudinal_waves.is_none() {
            return Err("axial bearings require explicit longitudinal-waves".into());
        }
        if let Some((resolution, config)) = &self.longitudinal_waves {
            resolution.validate()?;
            config.validate()?;
            if self.geometric_stretching.is_some() {
                return Err("longitudinal waves and geometric stretching cannot both supply strain energy".into());
            }
            modes = resolution.axial_modes.len().checked_mul(self.strings.len())
                .and_then(|axial| modes.checked_add(axial)).ok_or("complete longitudinal mode count overflow")?;
        }
        if let Some(budget) = self.prepared_linear {
            budget.validate()?;
            if modes > budget.max_modes {
                return Err("complete instrument exceeds the explicit prepared mode budget".into());
            }
        }
        if modes == 0 || modes > MAX_TOTAL_MODES {
            return Err(format!("complete instrument, including gas memory, exceeds {MAX_TOTAL_MODES} modes").into());
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
            "measured" => {
                let count = row.count(8)?;
                let mut axes = Vec::with_capacity(count);
                for _ in 0..count {
                    let axis = row.label()?;
                    let value = row.scalar()?;
                    if axes.iter().any(|(prior, _)| prior == &axis) {
                        return Err("measured string coordinates must have distinct axes".into());
                    }
                    axes.push((axis, value));
                }
                let properties = MeasuredStringProperties {
                    length: row.label()?, diameter: row.label()?, tension: row.label()?,
                    linear_density: row.label()?, bending_stiffness: row.label()?,
                };
                properties.requirements()?;
                if row.word()? != "no-intrinsic-loss" {
                    return Err("measured strings require explicit no-intrinsic-loss; assembly measurements do not identify a damping law".into());
                }
                row.end()?;
                Ok(WireSection::Measured { axes, properties })
            }
            "close-wound" => {
                let core_radius_m = row.positive()?;
                let pack = row.path()?;
                let wire_radius_m = row.positive()?;
                let pitch_m = row.positive()?;
                row.end()?;
                if ((core_radius_m + wire_radius_m) + wire_radius_m - outer_radius).abs()
                    > 8.0 * f64::EPSILON * outer_radius {
                    return Err("finished wound radius must equal core radius plus two wire radii".into());
                }
                Ok(WireSection::CloseWound { core_radius_m, pack, wire_radius_m, pitch_m })
            }
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
            _ => Err("section must be homogeneous, measured AXES AXIS VALUE ... LENGTH DIAMETER TENSION LINE_DENSITY EI no-intrinsic-loss, bonded CORE_RADIUS SLEEVE_COUNT, or close-wound CORE_RADIUS WRAP_PACK WIRE_RADIUS PITCH".into()),
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
    fn label(&mut self) -> Result<String> {
        let value = self.word()?;
        if value.len() > 256 || value.chars().any(char::is_control) {
            return Err("measurement axes and property names must be bounded labels".into());
        }
        Ok(value.to_owned())
    }
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

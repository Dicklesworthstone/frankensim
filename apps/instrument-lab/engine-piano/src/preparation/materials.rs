//! Canonical material packs reach the original live string/bridge/hammer owner.
//! Geometry is explicit; neither a material name nor a density implies radius.
//! Each string is a solid core of perfectly bonded concentric layers with
//! sourced bending loss (a Kelvin–Voigt viscosity or a loss tangent, each
//! layer its own), optionally under helical wraps that add only their wire
//! mass: the core carries the tension, EA, EI and loss.
//!
//! The construction also declares how tension is set and which sourced tensile
//! limit it is checked against. At a fixed pitch the stress `4 L^2 f^2 rho`
//! does not depend on the wire diameter, so a dense, weak metal cannot reach a
//! concert scale: every overstressed key refuses together, naming its required
//! stress, its limit and the highest pitch its section can reach at that
//! length. A pack without the limit claim refuses as missing data. Wrap
//! stiffness, interlayer slip and evolving temperature need their own physical
//! owners and are not inferred.

use std::{collections::{BTreeMap, BTreeSet}, sync::Arc};
use fs_matdb::{NormalizedMaterialCardPack, QueryPoint};
use fs_material::{section::{ConcentricElasticSection, ConcentricKelvinVoigtSection,
    ConcentricMaterialLayer, MAX_CONCENTRIC_LAYERS, loss_angle::ConcentricLossAngleSection,
    ted::ConcentricThermoelasticBending,
    strength::AxialTensileCapacity,
    wound::{HelicalWrap, HelicalWrapping, MAX_HELICAL_WRAPS, close_wound_pitch_m}},
    state_point::{DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY, MaterialPropertySelection,
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
    MaterialStatePointError, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point}};
use fs_qty::{Density, Dims, DynViscosity, Pressure};
use crate::refmods::geometry::{CoreLoss, Course};

pub const MAX_PACK_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_PACKS: usize = 128;
pub const HEADER: &str = "frankensim-piano-materials-v2";
pub const COLUMNS: &str = "key,layer,outer_radius_m,wire_diameter_m,pitch_m,material_index";
/// The superseded bonded-only format. It names no tensile limit, so the
/// mandatory yield check cannot run and it refuses as missing data.
pub const HEADER_V1: &str = "frankensim-piano-materials-v1";

/// Explicit browser operating point in absolute kelvin and cyclic hertz.
/// Axis names select original source coordinates; units and conventions are
/// fixed by this entry and cannot be inferred from an axis spelling.
pub fn operating_point(temperature_axis: &str, temperature_k: f64,
    frequency_axis: &str, frequency_hz: f64) -> Result<QueryPoint, String> {
    use fs_qty::{QuantitySpec, semantic::{FrequencyConvention, QuantityKind,
        SemanticType, ValueForm}};
    if [temperature_axis, frequency_axis].iter().any(|axis|
        axis.trim().is_empty() || axis.trim() != *axis || axis.len() > 128)
        || temperature_axis == frequency_axis
        || !temperature_k.is_finite() || temperature_k <= 0.0
        || !frequency_hz.is_finite() || frequency_hz <= 0.0 {
        return Err("declare distinct temperature/frequency axes and positive finite K/Hz values".into());
    }
    let typed = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    QueryPoint::new().with_quantity(temperature_axis, typed(QuantityKind::AbsoluteTemperature), temperature_k)
        .and_then(|point| point.with_quantity(frequency_axis,
            typed(QuantityKind::Frequency(FrequencyConvention::Cyclic)), frequency_hz))
        .map_err(|e| e.to_string())
}

/// How each course's tension is set once its material section is bound.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TensionPolicy {
    /// Re-solve tension so the new section's pinned stiff-string fundamental
    /// equals the installed course's fundamental (`target-pitch`).
    TargetPitch,
    /// Keep the installed tension; the pitch follows the new section
    /// (`installed-tension`).
    InstalledTension,
    /// Tension is this fraction, in (0, 1], of the core's axial capacity under
    /// the declared limit property (`yield-fraction,<x>`).
    YieldFraction(f64),
}
impl TensionPolicy {
    fn json(self) -> String {
        match self {
            Self::TargetPitch => "\"target-pitch\"".into(),
            Self::InstalledTension => "\"installed-tension\"".into(),
            Self::YieldFraction(x) => format!("{{\"yieldFraction\":{x}}}"),
        }
    }
}

/// Constitutive law of the cores' sourced bending loss, named by the
/// construction's optional `loss,` line; the bending-loss property passed to
/// the constructor is read under this law.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BendingLossLaw {
    /// Solid bending viscosity [Pa s] (`loss,kelvin-voigt`; also when no line).
    KelvinVoigt,
    /// Dimensionless loss tangent Im(E)/Re(E) (`loss,loss-angle`).
    LossAngle,
    /// Zener thermoelastic loss computed from each homogeneous core's E,
    /// density, expansion, heat capacity, conductivity and temperature, plus
    /// the property as a measured anelastic excess loss tangent inside its
    /// band (`loss,thermoelastic`). Without an excess claim, or outside its
    /// band, the excess is unknown and omitted: that decay is a lower bound.
    Thermoelastic,
}
impl BendingLossLaw {
    pub fn name(self) -> &'static str {
        match self {
            Self::KelvinVoigt => "kelvin-voigt",
            Self::LossAngle => "loss-angle",
            Self::Thermoelastic => "thermoelastic",
        }
    }
}

/// One course's tensile state: the most-tensioned unison member, the stress
/// in the core layer that reaches its limit first, and that limit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KeyStress {
    pub key: u8,
    /// Core layer index, innermost first, that governs the capacity.
    pub governing_layer: usize,
    pub max_member_tension_n: f64,
    pub required_stress_pa: f64,
    pub limit_pa: f64,
    /// Tension over capacity; above one the governing layer is overstressed.
    pub utilization: f64,
    pub capacity_n: f64,
    /// Highest pinned stiff-string fundamental this section reaches at the
    /// course's speaking length, at the capacity tension. For one material
    /// and negligible EI it is `sqrt(limit / rho) / (2 L)`.
    pub f_max_hz: f64,
}
impl KeyStress {
    fn json(&self) -> String {
        format!(concat!("{{\"key\":{},\"governingLayer\":{},\"maxMemberTensionN\":{},",
            "\"requiredStressPa\":{},\"limitPa\":{},\"utilization\":{},\"capacityN\":{},",
            "\"fMaxHz\":{}}}"), self.key, self.governing_layer, self.max_member_tension_n,
            self.required_stress_pa, self.limit_pa, self.utilization, self.capacity_n, self.f_max_hz)
    }
}

/// Every course's tensile state under one declared limit property.
#[derive(Clone, Debug, PartialEq)]
pub struct StrengthReport {
    pub limit_property: String,
    /// The labelled non-physical override was declared: overstressed strings
    /// play instead of refusing. Real strings at these stresses would break.
    pub non_physical_override: bool,
    pub keys: Vec<KeyStress>,
}
impl StrengthReport {
    pub fn overstressed(&self) -> impl Iterator<Item = &KeyStress> {
        self.keys.iter().filter(|key| key.utilization > 1.0)
    }

    /// The structured refusal: every overstressed key, with its numbers.
    pub fn refusal_json(&self) -> String {
        let keys = self.overstressed().map(KeyStress::json).collect::<Vec<_>>();
        format!("{{\"refusal\":\"overstressed-strings\",\"limitProperty\":\"{}\",\"keys\":{},\"overstressed\":[{}]}}",
            escape(&self.limit_property), self.keys.len(), keys.join(","))
    }
}

/// Admitted source owners and their actual section integrals. Cloning a live
/// engine retains these immutable packs, including original source receipts.
#[derive(Clone, Debug)]
pub struct StringMaterials {
    sections: BTreeMap<u8, ConcentricKelvinVoigtSection>,
    loss_angle_sections: BTreeMap<u8, ConcentricLossAngleSection>,
    /// Thermoelastic cores, with their measured excess when a source has one.
    thermoelastic_sections: BTreeMap<u8, (ConcentricThermoelasticBending, Option<ConcentricLossAngleSection>)>,
    /// Cores whose packs have no usable bending-loss claim: elastic only.
    unrated: BTreeMap<u8, ConcentricElasticSection>,
    /// The deferred missing-data refusal for `unrated`, raised only after the
    /// tensile check: a string that would break is reported first.
    missing_loss: Option<String>,
    loss: BendingLossLaw,
    wrappings: BTreeMap<u8, HelicalWrapping>,
    capacities: BTreeMap<u8, AxialTensileCapacity>,
    tension: TensionPolicy,
    limit_property: String,
    non_physical_override: bool,
    packs: Vec<Arc<NormalizedMaterialCardPack>>,
}

#[derive(Clone, Copy)]
enum Role { Core, Wrap }

struct Declaration {
    tension: TensionPolicy,
    loss: BendingLossLaw,
    limit_property: String,
    non_physical_override: bool,
    // Per key: core layers (outer radius, pack) then wraps (diameter, pitch or close, pack).
    cores: BTreeMap<u8, Vec<(f64, usize)>>,
    wraps: BTreeMap<u8, Vec<(f64, Option<f64>, usize)>>,
}

fn property_name(name: &str, what: &str) -> Result<String, String> {
    if name.trim().is_empty() || name.len() > 128 || name.trim() != name {
        return Err(format!("declare a nonempty trimmed {what} property of at most 128 bytes"));
    }
    Ok(name.into())
}

fn escape(text: &str) -> String {
    text.chars().map(|ch| if ch == '"' || ch == '\\' { '\'' } else if ch.is_control() { ' ' } else { ch }).collect()
}

fn parse_declaration(geometry: &str, packs: usize) -> Result<Declaration, String> {
    let mut lines = geometry.lines().map(|line| line.split('#').next().unwrap_or("").trim())
        .filter(|line| !line.is_empty());
    match lines.next() {
        Some(HEADER) => {}
        Some(HEADER_V1) => return Err(format!("{HEADER_V1} declares no tension policy or tensile-limit \
            property, so the mandatory yield check cannot run (missing data); use {HEADER}")),
        _ => return Err(format!("string material geometry requires {HEADER}")),
    }
    let tension = match lines.next().map(|line| line.split(',').map(str::trim).collect::<Vec<_>>()).as_deref() {
        Some(["tension", "target-pitch"]) => TensionPolicy::TargetPitch,
        Some(["tension", "installed-tension"]) => TensionPolicy::InstalledTension,
        Some(["tension", "yield-fraction", x]) => {
            let x = x.parse::<f64>().map_err(|_| "invalid yield fraction")?;
            if !x.is_finite() || x <= 0.0 || x > 1.0 { return Err("yield fraction must lie in (0, 1]".into()); }
            TensionPolicy::YieldFraction(x)
        }
        _ => return Err("second line must be tension,target-pitch | installed-tension | yield-fraction,<x>".into()),
    };
    let (limit_property, non_physical_override) =
        match lines.next().map(|line| line.split(',').map(str::trim).collect::<Vec<_>>()).as_deref() {
            Some(["strength", name]) => (property_name(name, "tensile-limit")?, false),
            Some(["strength", name, "non-physical-override"]) => (property_name(name, "tensile-limit")?, true),
            _ => return Err("third line must be strength,<tensile-limit property>[,non-physical-override]".into()),
        };
    let mut lines = lines.peekable();
    let loss = match lines.peek().map(|line| line.split(',').map(str::trim).collect::<Vec<_>>()).as_deref() {
        Some(["loss", "kelvin-voigt"]) => { lines.next(); BendingLossLaw::KelvinVoigt }
        Some(["loss", "loss-angle"]) => { lines.next(); BendingLossLaw::LossAngle }
        Some(["loss", "thermoelastic"]) => { lines.next(); BendingLossLaw::Thermoelastic }
        Some(["loss", ..]) => return Err("loss line must be loss,kelvin-voigt, loss,loss-angle or loss,thermoelastic".into()),
        _ => BendingLossLaw::KelvinVoigt,
    };
    if lines.next() != Some(COLUMNS) { return Err(format!("the line after strength (and loss) must be {COLUMNS}")); }
    let mut cores = BTreeMap::<u8, Vec<(f64, usize)>>::new();
    let mut wraps = BTreeMap::<u8, Vec<(f64, Option<f64>, usize)>>::new();
    for line in lines {
        let fields = line.split(',').map(str::trim).collect::<Vec<_>>();
        let [key, layer, radius, diameter, pitch, index] = fields.as_slice() else {
            return Err("string material row needs key,layer,outer_radius_m,wire_diameter_m,pitch_m,material_index".into());
        };
        let midi = key.parse::<u8>().map_err(|_| "invalid string material key")?;
        if !(21..=108).contains(&midi) { return Err("string material key must be 21..108".into()); }
        let index = index.parse::<usize>().map_err(|_| "invalid material pack index")?;
        if index >= packs { return Err("string material pack index is out of range".into()); }
        match *layer {
            "core" => {
                if !diameter.is_empty() || !pitch.is_empty() {
                    return Err(format!("key {midi}: a core row gives only outer_radius_m"));
                }
                if wraps.contains_key(&midi) { return Err(format!("key {midi}: core layers precede its wraps")); }
                let radius = radius.parse::<f64>().map_err(|_| "invalid string core outer radius")?;
                let course = cores.entry(midi).or_default();
                if course.len() >= MAX_CONCENTRIC_LAYERS {
                    return Err("string section exceeds the 128-layer input budget".into());
                }
                course.push((radius, index));
            }
            "wrap" => {
                if !radius.is_empty() { return Err(format!("key {midi}: a wrap row gives wire_diameter_m and pitch_m")); }
                if !cores.contains_key(&midi) { return Err(format!("key {midi}: a wrap needs a core beneath it")); }
                let diameter = diameter.parse::<f64>().map_err(|_| "invalid wrap wire diameter")?;
                let pitch = if *pitch == "close" { None }
                    else { Some(pitch.parse::<f64>().map_err(|_| "wrap pitch must be metres or close")?) };
                let course = wraps.entry(midi).or_default();
                if course.len() >= MAX_HELICAL_WRAPS {
                    return Err("string wrapping exceeds the 8-wrap input budget".into());
                }
                course.push((diameter, pitch, index));
            }
            _ => return Err(format!("key {midi}: layer must be core or wrap")),
        }
    }
    if cores.is_empty() { return Err("string material geometry has no sections".into()); }
    Ok(Declaration { tension, loss, limit_property, non_physical_override, cores, wraps })
}

impl StringMaterials {
    /// Read original normalized packs, then assign them by explicit zero-based
    /// index to each key's core layers (innermost first) and wraps. A core pack
    /// must resolve density, Young's modulus, the bending-loss property under
    /// the declared law and the declared tensile limit; a wrap pack must
    /// resolve density. Each pack
    /// resolves once, at one complete query, for every role it plays.
    /// Missing/ambiguous claims, wrong dimensions or unsupported state refuse.
    pub fn from_packs(bytes: &[&[u8]], geometry: &str, point: &QueryPoint,
        viscosity_property: &str) -> Result<Self, String> {
        if bytes.is_empty() || bytes.len() > MAX_PACKS {
            return Err("string materials require 1..128 canonical packs".into());
        }
        let mut size = 0usize;
        for pack in bytes {
            size = size.checked_add(pack.len()).ok_or("material pack byte overflow")?;
            if pack.is_empty() || size > MAX_PACK_BYTES {
                return Err("nonempty material packs must fit the 16 MiB input budget".into());
            }
        }
        if geometry.len() > super::MAX_SCALE_BYTES {
            return Err("string material geometry exceeds the 256 KiB input budget".into());
        }
        if viscosity_property.trim().is_empty() || viscosity_property.len() > 128
            || viscosity_property.trim() != viscosity_property
            || [DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY].contains(&viscosity_property) {
            return Err("declare a distinct nonempty solid bending-viscosity property".into());
        }
        let declaration = parse_declaration(geometry, bytes.len())?;
        let limit = declaration.limit_property.as_str();
        if [DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY, viscosity_property].contains(&limit) {
            return Err("the tensile-limit property must differ from density, modulus and viscosity".into());
        }
        let mut roles = vec![Vec::new(); bytes.len()];
        for layers in declaration.cores.values() {
            for &(_, index) in layers { roles[index].push(Role::Core); }
        }
        for layers in declaration.wraps.values() {
            for &(_, _, index) in layers { roles[index].push(Role::Wrap); }
        }
        let mut packs = Vec::with_capacity(bytes.len());
        let mut states = Vec::with_capacity(bytes.len());
        let mut loss_rated = Vec::with_capacity(bytes.len());
        let mut excess_measured = Vec::with_capacity(bytes.len());
        for (index, (bytes, roles)) in bytes.iter().zip(&roles).enumerate() {
            let pack = Arc::new(NormalizedMaterialCardPack::from_bytes(bytes)
                .map_err(|e| e.to_string())?);
            let core = roles.iter().any(|role| matches!(role, Role::Core));
            if roles.is_empty() { return Err(format!("string material pack {index} is not used by any layer")); }
            let mut requirements = vec![ScalarPropertyRequirement::try_new(DENSITY_PROPERTY, Density::DIMS,
                ScalarAdmissibility::StrictlyPositive)];
            if core {
                requirements.extend([
                    ScalarPropertyRequirement::try_new(YOUNG_MODULUS_PROPERTY, Pressure::DIMS,
                        ScalarAdmissibility::StrictlyPositive),
                    ScalarPropertyRequirement::try_new(viscosity_property, match declaration.loss {
                        BendingLossLaw::KelvinVoigt => DynViscosity::DIMS,
                        BendingLossLaw::LossAngle | BendingLossLaw::Thermoelastic => Dims::NONE,
                    }, ScalarAdmissibility::Finite),
                    ScalarPropertyRequirement::try_new(limit, Pressure::DIMS,
                        ScalarAdmissibility::StrictlyPositive),
                ]);
                if declaration.loss == BendingLossLaw::Thermoelastic {
                    requirements.extend([
                        ScalarPropertyRequirement::try_new(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
                            INVERSE_TEMPERATURE_DIMS, ScalarAdmissibility::Finite),
                        ScalarPropertyRequirement::try_new(SPECIFIC_HEAT_CAPACITY_PROPERTY,
                            SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility::StrictlyPositive),
                        ScalarPropertyRequirement::try_new(THERMAL_CONDUCTIVITY_PROPERTY,
                            THERMAL_CONDUCTIVITY_DIMS, ScalarAdmissibility::StrictlyPositive),
                    ]);
                }
            }
            let mut requirements = requirements.into_iter().collect::<Result<Vec<_>, _>>()
                .map_err(|e| e.to_string())?;
            let refuse = |e: MaterialStatePointError| match e {
                MaterialStatePointError::Query { ref property, .. } if property == limit => format!(
                    "string material pack {index} has no usable '{limit}' claim at this state: the \
                    mandatory yield check needs a sourced tensile limit (missing data): {e}"),
                e => format!("string material pack {index}: {e}"),
            };
            // Without loss data the core still has mass, stiffness and strength:
            // missing loss claims are dropped here and refuse only after the
            // tensile check. Under the thermoelastic law an absent excess claim
            // is not missing data: the excess is then unknown.
            let thermal = [LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, SPECIFIC_HEAT_CAPACITY_PROPERTY,
                THERMAL_CONDUCTIVITY_PROPERTY];
            let mut dropped = Vec::new();
            let state = loop {
                match resolve_material_state_point(pack.card(), point, &requirements,
                    MaterialPropertySelection::SingleClaimOnly) {
                    Err(MaterialStatePointError::Query { property, .. })
                        if core && (property == viscosity_property || (declaration.loss == BendingLossLaw::Thermoelastic
                            && thermal.contains(&property.as_str()))) => {
                        requirements.retain(|requirement| requirement.name() != property);
                        dropped.push(property);
                    }
                    state => break state.map_err(refuse)?,
                }
            };
            let rated = match declaration.loss {
                BendingLossLaw::Thermoelastic => !dropped.iter().any(|property| thermal.contains(&property.as_str())),
                _ => dropped.is_empty(),
            };
            packs.push(pack);
            states.push(state);
            loss_rated.push(rated);
            excess_measured.push(!dropped.iter().any(|property| property == viscosity_property));
        }
        Self::from_states(declaration, &states, &loss_rated, &excess_measured, viscosity_property, packs)
    }

    fn from_states(declaration: Declaration, states: &[ResolvedMaterialStatePoint], loss_rated: &[bool],
        excess_measured: &[bool], viscosity_property: &str, packs: Vec<Arc<NormalizedMaterialCardPack>>) -> Result<Self, String> {
        let Declaration { tension, loss, limit_property, non_physical_override, cores, mut wraps } = declaration;
        let mut sections = BTreeMap::new();
        let mut loss_angle_sections = BTreeMap::new();
        let mut thermoelastic_sections = BTreeMap::new();
        let mut unrated = BTreeMap::new();
        let mut wrappings = BTreeMap::new();
        let mut capacities = BTreeMap::new();
        for (key, layers) in cores {
            let failed = |e: MaterialStatePointError| format!("string material key {key}: {e}");
            let rated = layers.iter().all(|&(_, index)| loss_rated[index]);
            let measured = layers.iter().all(|&(_, index)| excess_measured[index]);
            let layers = layers.into_iter().map(|(outer_radius_m, index)|
                ConcentricMaterialLayer { outer_radius_m, material: states[index].clone() }).collect();
            let elastic = ConcentricElasticSection::try_new(0.0, layers).map_err(failed)?;
            capacities.insert(key, elastic.axial_tensile_capacity(&limit_property).map_err(failed)?);
            if let Some(rows) = wraps.remove(&key) {
                let mut under = 2.0 * elastic.outer_radius_m();
                let mut helical = Vec::with_capacity(rows.len());
                for (wire_diameter_m, pitch, index) in rows {
                    let pitch_m = match pitch {
                        Some(pitch) => pitch,
                        None => close_wound_pitch_m(under, wire_diameter_m).map_err(failed)?,
                    };
                    helical.push(HelicalWrap { wire_diameter_m, pitch_m, material: states[index].clone() });
                    under += 2.0 * wire_diameter_m;
                }
                wrappings.insert(key, HelicalWrapping::try_new(&elastic, helical).map_err(failed)?);
            }
            if !rated {
                unrated.insert(key, elastic);
                continue;
            }
            match loss {
                BendingLossLaw::KelvinVoigt => {
                    sections.insert(key, elastic.with_kelvin_voigt_bending(viscosity_property).map_err(failed)?);
                }
                BendingLossLaw::LossAngle => {
                    loss_angle_sections.insert(key, elastic.with_loss_angle_bending(viscosity_property).map_err(failed)?);
                }
                BendingLossLaw::Thermoelastic => {
                    let excess = if measured {
                        Some(elastic.clone().with_loss_angle_bending(viscosity_property).map_err(failed)?)
                    } else { None };
                    thermoelastic_sections.insert(key, (elastic.with_thermoelastic_bending().map_err(failed)?, excess));
                }
            }
        }
        let missing_loss = (!unrated.is_empty()).then(|| {
            let keys = unrated.keys().map(u8::to_string).collect::<Vec<_>>().join(",");
            let unrated_packs = loss_rated.iter().enumerate().filter(|(_, rated)| !**rated)
                .map(|(index, _)| index.to_string()).collect::<Vec<_>>().join(",");
            let data = match loss {
                BendingLossLaw::Thermoelastic => "thermal expansion, heat capacity and conductivity".to_string(),
                law => format!("'{viscosity_property}' {}", law.name()),
            };
            format!("string material packs [{unrated_packs}] have no usable {data} \
                claims at this state, so the cores of keys [{keys}] have no sourced bending loss (missing data)")
        });
        Ok(Self { sections, loss_angle_sections, thermoelastic_sections, unrated, missing_loss, loss, wrappings,
            capacities, tension, limit_property, non_physical_override, packs })
    }

    /// Bounded WASM framing of the unchanged canonical bytes; this is not an
    /// alternate material serialization. Every byte must belong to one pack.
    pub fn from_packed_bytes(bytes: &[u8], lengths: &[u32], geometry: &str,
        point: &QueryPoint, viscosity_property: &str) -> Result<Self, String> {
        if bytes.len() > MAX_PACK_BYTES || lengths.is_empty() || lengths.len() > MAX_PACKS {
            return Err("string material packet exceeds the pack count or byte budget".into());
        }
        let mut offset = 0usize;
        let mut packs = Vec::with_capacity(lengths.len());
        for &length in lengths {
            let end = offset.checked_add(length as usize).ok_or("material framing overflow")?;
            packs.push(bytes.get(offset..end).ok_or("truncated material pack packet")?);
            offset = end;
        }
        if offset != bytes.len() { return Err("trailing bytes after string material packs".into()); }
        Self::from_packs(&packs, geometry, point, viscosity_property)
    }

    /// Bonded Kelvin–Voigt core sections, which carry tension, EA, EI and
    /// bending loss (empty under the loss-angle law).
    pub fn sections(&self) -> &BTreeMap<u8, ConcentricKelvinVoigtSection> { &self.sections }
    /// Bonded loss-angle core sections (empty under the Kelvin–Voigt law).
    pub fn loss_angle_sections(&self) -> &BTreeMap<u8, ConcentricLossAngleSection> { &self.loss_angle_sections }
    pub fn bending_loss_law(&self) -> BendingLossLaw { self.loss }
    /// Thermoelastic cores and their measured excess, if any (empty unless
    /// the thermoelastic law was declared).
    pub fn thermoelastic_sections(&self)
        -> &BTreeMap<u8, (ConcentricThermoelasticBending, Option<ConcentricLossAngleSection>)> {
        &self.thermoelastic_sections
    }
    /// Helical wraps of the wound keys; they add only line density.
    pub fn wrappings(&self) -> &BTreeMap<u8, HelicalWrapping> { &self.wrappings }
    /// Each core's axial capacity under the declared tensile limit.
    pub fn capacities(&self) -> &BTreeMap<u8, AxialTensileCapacity> { &self.capacities }
    pub fn tension_policy(&self) -> TensionPolicy { self.tension }
    pub fn limit_property(&self) -> &str { &self.limit_property }
    pub fn non_physical_override(&self) -> bool { self.non_physical_override }
    pub fn packs(&self) -> &[Arc<NormalizedMaterialCardPack>] { &self.packs }

    /// The deferred refusal when some cores have no sourced bending loss.
    pub fn missing_loss(&self) -> Option<&str> { self.missing_loss.as_deref() }

    /// Finished outer radius of a key's string: its outermost wrap, else its core.
    pub fn outer_radius_m(&self, key: u8) -> Option<f64> {
        self.wrappings.get(&key).map(HelicalWrapping::outer_radius_m)
            .or_else(|| self.elastic(key).map(ConcentricElasticSection::outer_radius_m))
    }

    fn elastic(&self, key: u8) -> Option<&ConcentricElasticSection> {
        self.core(key).map(CoreLoss::elastic).or_else(|| self.unrated.get(&key))
    }

    fn core(&self, key: u8) -> Option<CoreLoss<'_>> {
        self.sections.get(&key).map(CoreLoss::KelvinVoigt)
            .or_else(|| self.loss_angle_sections.get(&key).map(CoreLoss::LossAngle))
            .or_else(|| self.thermoelastic_sections.get(&key)
                .map(|(section, excess)| CoreLoss::Thermoelastic(section, excess.as_ref())))
    }

    fn covers(&self, courses: &[Course]) -> Result<(), String> {
        let keys = courses.iter().map(|course| course.midi).collect::<BTreeSet<_>>();
        if keys.len() != courses.len() || keys != self.capacities.keys().copied().collect() {
            return Err("string materials must cover exactly every selected scale key".into());
        }
        Ok(())
    }

    /// Material coefficients only; the course's own tension is kept.
    fn bind(&self, course: &Course) -> Result<Course, String> {
        let Some(core) = self.core(course.midi) else {
            return Err(self.missing_loss.clone().unwrap_or_else(|| "string material key has no core section".into()));
        };
        course.with_loss_section(core, self.wrappings.get(&course.midi))
            .map_err(|e| format!("string material key {}: {e}", course.midi))
    }

    /// Bind the complete selected scale, set every tension by the declared
    /// policy, then check every string against its sourced tensile limit.
    /// Length, hammer, unison and duplex geometry survive; μ, EA, EI and ηI
    /// use source integrals. Overstressed strings refuse together unless the
    /// labelled non-physical override was declared; only then does missing
    /// bending-loss data refuse.
    pub fn apply(&self, courses: &[Course]) -> Result<Vec<Course>, String> {
        self.covers(courses)?;
        let tensioned = courses.iter().map(|course| self.tensioned(course)).collect::<Result<Vec<_>, String>>()?;
        self.check_strength(&tensioned)?;
        if let Some(missing) = &self.missing_loss { return Err(missing.clone()); }
        courses.iter().zip(&tensioned).map(|(course, tensioned)| {
            let mut changed = self.bind(course)?;
            changed.tension_n = tensioned.tension_n;
            changed.validate().map_err(|e| format!("string material key {}: {e}", course.midi))?;
            Ok(changed)
        }).collect()
    }

    /// The course with its section's mass and stiffness and its policy
    /// tension, before any loss binding. Same arithmetic as the bound course.
    fn tensioned(&self, course: &Course) -> Result<Course, String> {
        let elastic = self.elastic(course.midi).ok_or("string material key has no core section")?;
        let mut changed = Course {
            linear_density_kg_m: match self.wrappings.get(&course.midi) {
                Some(wrapping) => elastic.line_density_kg_m() + wrapping.line_density_kg_m(),
                None => elastic.line_density_kg_m(),
            },
            flexural_rigidity_nm2: elastic.bending_stiffness_n_m2(),
            bonded_section: None,
            ..*course
        };
        {
            changed.tension_n = match self.tension {
                TensionPolicy::InstalledTension => course.tension_n,
                TensionPolicy::TargetPitch => {
                    // Invert the pinned stiff-string law for the installed
                    // fundamental: T = mu (2 L f)^2 - EI (pi / L)^2.
                    let target = course.partial_hz(1, course.tension_n);
                    let k = std::f64::consts::PI / changed.length_m;
                    changed.linear_density_kg_m * (2.0 * changed.length_m * target).powi(2)
                        - changed.flexural_rigidity_nm2 * k * k
                }
                TensionPolicy::YieldFraction(x) => x * self.capacities[&course.midi].capacity_n(),
            };
        }
        changed.validate().map_err(|e| format!("string material key {}: {e}", course.midi))?;
        Ok(changed)
    }

    /// Tensile state of every course at its actual tensions, including the
    /// detuned unison members.
    pub fn strength_report(&self, courses: &[Course]) -> Result<StrengthReport, String> {
        self.covers(courses)?;
        let mut keys = Vec::with_capacity(courses.len());
        for course in courses {
            let capacity = &self.capacities[&course.midi];
            let mut tension = 0.0_f64;
            for member in 0..course.unison {
                let cents = (member as f64 - 0.5 * (course.unison - 1) as f64) * course.detune_cents;
                tension = tension.max(course.tension_at_cents(cents)?);
            }
            let layer = capacity.governing_layer();
            let fail = |e: MaterialStatePointError| format!("string material key {}: {e}", course.midi);
            let f_max_hz = course.partial_hz(1, capacity.capacity_n());
            if !f_max_hz.is_finite() { return Err(format!("key {}: unrepresentable highest pitch", course.midi)); }
            keys.push(KeyStress {
                key: course.midi, governing_layer: layer, max_member_tension_n: tension,
                required_stress_pa: capacity.layer_stress_pa(layer, tension).map_err(fail)?,
                limit_pa: capacity.layer_limits_pa()[layer],
                utilization: capacity.utilization(tension).map_err(fail)?,
                capacity_n: capacity.capacity_n(), f_max_hz,
            });
        }
        Ok(StrengthReport { limit_property: self.limit_property.clone(),
            non_physical_override: self.non_physical_override, keys })
    }

    fn check_strength(&self, courses: &[Course]) -> Result<(), String> {
        let report = self.strength_report(courses)?;
        if !self.non_physical_override && report.overstressed().next().is_some() {
            return Err(format!("strings exceed their sourced tensile limit: {}", report.refusal_json()));
        }
        Ok(())
    }

    /// Admit courses built elsewhere (for example retuned after `apply`):
    /// their coefficients must equal their bound sections, and their actual
    /// tensions must pass the same tensile check.
    pub(crate) fn verify_courses(&self, courses: &[Course]) -> Result<(), String> {
        self.covers(courses)?;
        for course in courses {
            if self.bind(course)? != *course {
                return Err("live string coefficients differ from their bound material sections".into());
            }
        }
        self.check_strength(courses)
    }

    /// Actual source/state identities and coefficients installed in playback.
    pub fn json(&self) -> String {
        let sections = self.capacities.keys().map(|key| {
            let elastic = self.elastic(*key).expect("every capacity has its core section");
            // A core without sourced loss data reports null loss, never zero.
            let (viscous, hysteretic, band) = match self.core(*key) {
                Some(CoreLoss::KelvinVoigt(section)) => (section.viscous_stiffness_n_m2_s().to_string(),
                    "0".to_string(), section.omega_band_rad_s()),
                Some(CoreLoss::LossAngle(section)) => ("0".to_string(),
                    section.hysteretic_stiffness_n_m2().to_string(), section.omega_band_rad_s()),
                Some(CoreLoss::Thermoelastic(_, Some(excess))) => ("0".to_string(),
                    excess.hysteretic_stiffness_n_m2().to_string(), excess.omega_band_rad_s()),
                Some(CoreLoss::Thermoelastic(_, None)) | None =>
                    ("null".to_string(), "null".to_string(), (f64::NAN, f64::NAN)),
            };
            // Computed Zener part and whether an anelastic excess was measured.
            let thermoelastic = match self.core(*key) {
                Some(CoreLoss::Thermoelastic(section, excess)) => format!(
                    "{{\"rigidityNm2\":{},\"baseTimeS\":{},\"temperatureK\":{},\"excessLoss\":\"{}\"}}",
                    section.hysteretic_scale_n_m2(), section.base_time_s(), section.law().t0,
                    if excess.is_some() { "measured-in-band" } else { "unknown" }),
                _ => "null".to_string(),
            };
            let band = if band.0.is_nan() { "null".to_string() } else { format!("[{},{}]", band.0, band.1) };
            let identities = elastic.layers().iter().map(|layer|
                format!("\"{}\"", layer.material.identity().to_hex()))
                .collect::<Vec<_>>().join(",");
            let radii = elastic.layers().iter().map(|layer| layer.outer_radius_m)
                .collect::<Vec<_>>();
            let capacity = &self.capacities[key];
            let wraps = self.wrappings.get(key).map_or_else(|| "[]".into(), |wrapping| {
                let rows = wrapping.wraps().iter().zip(wrapping.wrap_line_density_kg_m()).map(|(wrap, mass)|
                    format!("{{\"stateIdentity\":\"{}\",\"wireDiameterM\":{},\"pitchM\":{},\"lineDensityKgM\":{}}}",
                        wrap.material.identity().to_hex(), wrap.wire_diameter_m, wrap.pitch_m, mass))
                    .collect::<Vec<_>>();
                format!("[{}]", rows.join(","))
            });
            format!(concat!("{{\"key\":{},\"stateIdentities\":[{}],\"outerRadiiM\":{:?},",
                "\"lineDensityKgM\":{},\"axialRigidityN\":{},\"flexuralRigidityNm2\":{},",
                "\"viscousRigidityNm2S\":{},\"hystereticRigidityNm2\":{},\"omegaBandRadS\":{},\"thermoelastic\":{},\"wraps\":{},",
                "\"outerRadiusM\":{},\"tensileCapacityN\":{},\"governingLayer\":{},\"layerLimitsPa\":{:?}}}"),
                key, identities, radii, elastic.line_density_kg_m()
                    + self.wrappings.get(key).map_or(0.0, HelicalWrapping::line_density_kg_m),
                elastic.axial_stiffness_n(), elastic.bending_stiffness_n_m2(),
                viscous, hysteretic, band, thermoelastic, wraps,
                self.outer_radius_m(*key).unwrap_or(0.0), capacity.capacity_n(),
                capacity.governing_layer(), capacity.layer_limits_pa())
        }).collect::<Vec<_>>().join(",");
        let receipt = match self.loss {
            BendingLossLaw::Thermoelastic => "thermoelastic loss is computed at every mode; a measured anelastic \
                excess is added only inside its source band; elsewhere intrinsic anelastic loss is omitted, so the \
                decay there is a lower bound",
            _ => "measured bending loss only; a mode outside its source band refuses",
        };
        format!(concat!("{{\"model\":\"bonded-core-kelvin-voigt-helical-wraps\",\"importedPacks\":{},",
            "\"tensionPolicy\":{},\"bendingLoss\":\"{}\",\"lossReceipt\":\"{}\",\"missingLossData\":{},",
            "\"limitProperty\":\"{}\",\"nonPhysicalOverride\":{},\"sections\":[{}]}}"),
            self.packs.len(), self.tension.json(), self.loss.name(), receipt, self.missing_loss.is_some(),
            escape(&self.limit_property), self.non_physical_override, sections)
    }
}

//! Canonical material packs reach the original live string/bridge/hammer owner.
//! Geometry is explicit; neither a material name nor a density implies radius.
//! Each string is a solid core of perfectly bonded concentric layers with
//! sourced Kelvin–Voigt bending loss, optionally under helical wraps that add
//! only their wire mass: the core carries the tension, EA, EI and loss.
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
    ConcentricMaterialLayer, MAX_CONCENTRIC_LAYERS, strength::AxialTensileCapacity,
    wound::{HelicalWrap, HelicalWrapping, MAX_HELICAL_WRAPS, close_wound_pitch_m}},
    state_point::{DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY, MaterialPropertySelection,
    MaterialStatePointError, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point}};
use fs_qty::{Density, DynViscosity, Pressure};
use crate::refmods::geometry::Course;

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
    if lines.next() != Some(COLUMNS) { return Err(format!("fourth line must be {COLUMNS}")); }
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
    Ok(Declaration { tension, limit_property, non_physical_override, cores, wraps })
}

impl StringMaterials {
    /// Read original normalized packs, then assign them by explicit zero-based
    /// index to each key's core layers (innermost first) and wraps. A core pack
    /// must resolve density, Young's modulus, the bending viscosity and the
    /// declared tensile limit; a wrap pack must resolve density. Each pack
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
                    ScalarPropertyRequirement::try_new(viscosity_property, DynViscosity::DIMS,
                        ScalarAdmissibility::Finite),
                    ScalarPropertyRequirement::try_new(limit, Pressure::DIMS,
                        ScalarAdmissibility::StrictlyPositive),
                ]);
            }
            let requirements = requirements.into_iter().collect::<Result<Vec<_>, _>>()
                .map_err(|e| e.to_string())?;
            let state = resolve_material_state_point(pack.card(), point, &requirements,
                MaterialPropertySelection::SingleClaimOnly).map_err(|e| match e {
                    MaterialStatePointError::Query { ref property, .. } if property == limit => format!(
                        "string material pack {index} has no usable '{limit}' claim at this state: the \
                        mandatory yield check needs a sourced tensile limit (missing data): {e}"),
                    e => format!("string material pack {index}: {e}"),
                })?;
            packs.push(pack);
            states.push(state);
        }
        Self::from_states(declaration, &states, viscosity_property, packs)
    }

    fn from_states(declaration: Declaration, states: &[ResolvedMaterialStatePoint],
        viscosity_property: &str, packs: Vec<Arc<NormalizedMaterialCardPack>>) -> Result<Self, String> {
        let Declaration { tension, limit_property, non_physical_override, cores, mut wraps } = declaration;
        let mut sections = BTreeMap::new();
        let mut wrappings = BTreeMap::new();
        let mut capacities = BTreeMap::new();
        for (key, layers) in cores {
            let failed = |e: MaterialStatePointError| format!("string material key {key}: {e}");
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
            sections.insert(key, elastic.with_kelvin_voigt_bending(viscosity_property).map_err(failed)?);
        }
        Ok(Self { sections, wrappings, capacities, tension, limit_property, non_physical_override, packs })
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

    /// Bonded core sections, which carry tension, EA, EI and bending loss.
    pub fn sections(&self) -> &BTreeMap<u8, ConcentricKelvinVoigtSection> { &self.sections }
    /// Helical wraps of the wound keys; they add only line density.
    pub fn wrappings(&self) -> &BTreeMap<u8, HelicalWrapping> { &self.wrappings }
    /// Each core's axial capacity under the declared tensile limit.
    pub fn capacities(&self) -> &BTreeMap<u8, AxialTensileCapacity> { &self.capacities }
    pub fn tension_policy(&self) -> TensionPolicy { self.tension }
    pub fn limit_property(&self) -> &str { &self.limit_property }
    pub fn non_physical_override(&self) -> bool { self.non_physical_override }
    pub fn packs(&self) -> &[Arc<NormalizedMaterialCardPack>] { &self.packs }

    /// Finished outer radius of a key's string: its outermost wrap, else its core.
    pub fn outer_radius_m(&self, key: u8) -> Option<f64> {
        self.wrappings.get(&key).map(HelicalWrapping::outer_radius_m)
            .or_else(|| self.sections.get(&key).map(|section| section.elastic().outer_radius_m()))
    }

    fn covers(&self, courses: &[Course]) -> Result<(), String> {
        let keys = courses.iter().map(|course| course.midi).collect::<BTreeSet<_>>();
        if keys.len() != courses.len() || keys != self.sections.keys().copied().collect() {
            return Err("string materials must cover exactly every selected scale key".into());
        }
        Ok(())
    }

    /// Material coefficients only; the course's own tension is kept.
    fn bind(&self, course: &Course) -> Result<Course, String> {
        let core = &self.sections[&course.midi];
        match self.wrappings.get(&course.midi) {
            Some(wrapping) => course.with_wound_kelvin_voigt_section(core, wrapping),
            None => course.with_bonded_kelvin_voigt_section(core),
        }.map_err(|e| format!("string material key {}: {e}", course.midi))
    }

    /// Bind the complete selected scale, set every tension by the declared
    /// policy, then check every string against its sourced tensile limit.
    /// Length, hammer, unison and duplex geometry survive; μ, EA, EI and ηI
    /// use source integrals. Overstressed strings refuse together unless the
    /// labelled non-physical override was declared.
    pub fn apply(&self, courses: &[Course]) -> Result<Vec<Course>, String> {
        self.covers(courses)?;
        let bound = courses.iter().map(|course| {
            let mut changed = self.bind(course)?;
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
            changed.validate().map_err(|e| format!("string material key {}: {e}", course.midi))?;
            Ok(changed)
        }).collect::<Result<Vec<_>, String>>()?;
        self.check_strength(&bound)?;
        Ok(bound)
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
        let sections = self.sections.iter().map(|(key, section)| {
            let elastic = section.elastic();
            let identities = elastic.layers().iter().map(|layer|
                format!("\"{}\"", layer.material.identity().to_hex()))
                .collect::<Vec<_>>().join(",");
            let radii = elastic.layers().iter().map(|layer| layer.outer_radius_m)
                .collect::<Vec<_>>();
            let band = section.omega_band_rad_s();
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
                "\"viscousRigidityNm2S\":{},\"omegaBandRadS\":[{},{}],\"wraps\":{},",
                "\"outerRadiusM\":{},\"tensileCapacityN\":{},\"governingLayer\":{},\"layerLimitsPa\":{:?}}}"),
                key, identities, radii, elastic.line_density_kg_m()
                    + self.wrappings.get(key).map_or(0.0, HelicalWrapping::line_density_kg_m),
                elastic.axial_stiffness_n(), elastic.bending_stiffness_n_m2(),
                section.viscous_stiffness_n_m2_s(), band.0, band.1, wraps,
                self.outer_radius_m(*key).unwrap_or(0.0), capacity.capacity_n(),
                capacity.governing_layer(), capacity.layer_limits_pa())
        }).collect::<Vec<_>>().join(",");
        format!(concat!("{{\"model\":\"bonded-core-kelvin-voigt-helical-wraps\",\"importedPacks\":{},",
            "\"tensionPolicy\":{},\"limitProperty\":\"{}\",\"nonPhysicalOverride\":{},\"sections\":[{}]}}"),
            self.packs.len(), self.tension.json(), escape(&self.limit_property),
            self.non_physical_override, sections)
    }
}

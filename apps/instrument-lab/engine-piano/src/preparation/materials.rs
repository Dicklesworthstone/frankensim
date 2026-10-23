//! Canonical material packs reach the original live string/bridge/hammer owner.
//! Geometry is explicit; neither a material name nor a density implies radius.
//! Each string is a solid core of perfectly bonded concentric layers with
//! sourced bending loss (a Kelvin–Voigt viscosity or a loss tangent, each
//! layer its own), optionally under helical wraps that add only their wire
//! mass: the core carries the tension, EA, EI and loss.
//!
//! The construction also declares how tension is set and which sourced tensile
//! limit it is checked against. Drawing work-hardens wire, so a strength claim
//! holds at the wire diameters it states: each core's limit is evaluated at the
//! core's own drawn (outer) diameter, from its layer packs or from one named
//! strength pack. A core outside the stated diameters refuses as missing data
//! unless the construction declares `extrapolate-diameter` (the nearest stated
//! diameter's value is taken); a claim that states no diameter at all refuses
//! unless it declares `diameter-free`, and then it cannot set a playing tension
//! (`yield-fraction` refuses). Each receipt prints the resolved limit, the rows
//! it came from and what it is not. Core and wrap diameters inferred rather than
//! measured or cited are declared (`geometry,back-solved,<first>-<last>`), and
//! every receipt and strength report repeats that caveat. At a fixed pitch the stress `4 L^2 f^2 rho`
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
    MaterialStatePointError, ResolvedMaterialStatePoint, ResolvedScalarProperty, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point}};
use fs_matdb::{EvaluationDecision, PropertyValue};
use fs_qty::{Density, Dims, DynViscosity, Length, Pressure, QuantitySpec};
use crate::refmods::geometry::{CoreLoss, Course};

pub const MAX_PACK_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_PACKS: usize = 128;
pub const HEADER: &str = "frankensim-piano-materials-v2";
pub const COLUMNS: &str = "key,layer,outer_radius_m,wire_diameter_m,pitch_m,material_index";
/// Query axis of a strength claim's wire diameter [m] (dimensional length).
pub const WIRE_DIAMETER_AXIS: &str = "wire_diameter";
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
    /// Wire diameter [m] at which the governing limit was evaluated; `None`
    /// when its claim states no diameter (applied only under `diameter-free`).
    pub limit_evaluated_at_m: Option<f64>,
    /// The governing limit was extrapolated beyond its stated diameters.
    pub limit_extrapolated: bool,
}
impl KeyStress {
    fn json(&self) -> String {
        let at = self.limit_evaluated_at_m.map_or_else(|| "null".to_string(), |d| d.to_string());
        format!(concat!("{{\"key\":{},\"governingLayer\":{},\"maxMemberTensionN\":{},",
            "\"requiredStressPa\":{},\"limitPa\":{},\"utilization\":{},\"capacityN\":{},",
            "\"fMaxHz\":{},\"limitAtDiameterM\":{},\"limitExtrapolated\":{}}}"), self.key, self.governing_layer,
            self.max_member_tension_n, self.required_stress_pa, self.limit_pa, self.utilization, self.capacity_n,
            self.f_max_hz, at, self.limit_extrapolated)
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
    /// What the limits are not (extrapolated, or not diameter-resolved), once
    /// per claim, in key order.
    pub limit_caveats: Vec<String>,
    /// The declared back-solved geometry of some of these keys, if any.
    pub geometry_caveat: Option<String>,
}
impl StrengthReport {
    pub fn overstressed(&self) -> impl Iterator<Item = &KeyStress> {
        self.keys.iter().filter(|key| key.utilization > 1.0)
    }

    /// The structured refusal: every overstressed key, with its numbers.
    pub fn refusal_json(&self) -> String {
        let keys = self.overstressed().map(KeyStress::json).collect::<Vec<_>>();
        let caveats = self.limit_caveats.iter().map(|c| format!("\"{}\"", escape(c))).collect::<Vec<_>>();
        let geometry = self.geometry_caveat.as_ref().map_or_else(|| "null".to_string(), |c| format!("\"{}\"", escape(c)));
        format!(concat!("{{\"refusal\":\"overstressed-strings\",\"limitProperty\":\"{}\",\"limitCaveats\":[{}],",
            "\"geometryCaveat\":{},\"keys\":{},\"overstressed\":[{}]}}"),
            escape(&self.limit_property), caveats.join(","), geometry, self.keys.len(), keys.join(","))
    }
}

/// Where one core layer's tensile limit came from.
#[derive(Clone, Debug, PartialEq)]
pub struct LimitReceipt {
    /// Pack whose claim supplied the limit.
    pub pack: usize,
    /// The drawn wire's diameter [m]: the core's outer diameter.
    pub diameter_m: f64,
    /// Diameter [m] at which the claim was evaluated: the core's own unless
    /// extrapolated to the nearest stated diameter; `None` when the claims
    /// state no diameter.
    pub evaluated_at_m: Option<f64>,
    /// The claims' stated wire-diameter range [m], if any.
    pub domain_m: Option<(f64, f64)>,
    /// The resolved limit [Pa].
    pub limit_pa: f64,
    /// How the claim was evaluated (tabulated row, linear rule, constant).
    pub decision: String,
    /// The selected claim's rows the limit came from, as (abscissa SI, limit
    /// Pa): both bracketing rows of a linear rule, the one row hit exactly,
    /// none for a scalar.
    pub rows: Vec<(f64, f64)>,
    /// The core lies outside the stated diameters or none is stated: the
    /// limit is not the claim's value at this core's diameter.
    pub extrapolated: bool,
    /// What the limit is not, when `extrapolated`.
    pub caveat: Option<String>,
}
impl LimitReceipt {
    fn json(&self) -> String {
        let opt = |v: Option<f64>| v.map_or_else(|| "null".to_string(), |v| v.to_string());
        let domain = self.domain_m.map_or_else(|| "null".to_string(), |(lo, hi)| format!("[{lo},{hi}]"));
        let rows = self.rows.iter().map(|(x, y)| format!("[{x},{y}]")).collect::<Vec<_>>().join(",");
        let caveat = self.caveat.as_ref().map_or_else(|| "null".to_string(), |c| format!("\"{}\"", escape(c)));
        format!(concat!("{{\"pack\":{},\"diameterM\":{},\"evaluatedAtM\":{},\"domainM\":{},\"limitPa\":{},",
            "\"decision\":\"{}\",\"rows\":[{}],\"extrapolated\":{},\"caveat\":{}}}"),
            self.pack, self.diameter_m, opt(self.evaluated_at_m), domain, self.limit_pa, escape(&self.decision),
            rows, self.extrapolated, caveat)
    }
}

/// The two strength flags a construction may declare for claims that do not
/// state this core's diameter.
#[derive(Clone, Copy)]
struct DiameterFlags { extrapolate: bool, diameter_free: bool }

/// Evaluate a tensile-limit claim of `pack` at a core's drawn diameter.
fn resolve_limit(pack: &NormalizedMaterialCardPack, point: &QueryPoint, limit: &str, index: usize, key: u8,
    diameter_m: f64, flags: DiameterFlags) -> Result<(ResolvedScalarProperty, LimitReceipt), String> {
    let claims = pack.card().claims().claims_for(limit);
    if claims.is_empty() {
        return Err(format!("string material pack {index} has no usable '{limit}' claim at this state: the \
            mandatory yield check needs a sourced tensile limit (missing data)"));
    }
    let domain = claims.iter().filter_map(|(_, claim)| claim.validity.bounds().get(WIRE_DIAMETER_AXIS).copied())
        .reduce(|(a, b), (c, d)| (a.min(c), b.max(d)));
    let mm = |m: f64| m * 1e3;
    let query_at = match domain {
        None if !flags.diameter_free => return Err(format!("string material key {key}: pack {index}'s '{limit}' \
            claims state no wire diameter, so none is a strength of this {:.4} mm core (missing data); declare \
            strength,{limit},diameter-free to apply them at every diameter, printed in the receipt", mm(diameter_m))),
        None => None,
        Some((lo, hi)) if (lo..=hi).contains(&diameter_m) => Some(diameter_m),
        Some((lo, hi)) if !flags.extrapolate => return Err(format!("string material key {key}: its {:.4} mm core lies \
            outside the {:.4}-{:.4} mm wire diameters of pack {index}'s '{limit}' claims, so its strength is \
            unknown (missing data); declare strength,{limit},extrapolate-diameter to take the nearest stated \
            diameter's value, printed in the receipt", mm(diameter_m), mm(lo), mm(hi))),
        Some((lo, hi)) => Some(diameter_m.clamp(lo, hi)),
    };
    let at = match query_at {
        Some(d) => point.clone().with_quantity(WIRE_DIAMETER_AXIS, QuantitySpec::dimensional(Length::DIMS), d)
            .map_err(|e| format!("string material key {key}: {e}"))?,
        None => point.clone(),
    };
    let requirement = ScalarPropertyRequirement::try_new(limit, Pressure::DIMS, ScalarAdmissibility::StrictlyPositive)
        .map_err(|e| e.to_string())?;
    let state = resolve_material_state_point(pack.card(), &at, &[requirement], MaterialPropertySelection::SingleClaimOnly)
        .map_err(|e| match query_at {
            Some(d) => format!("string material key {key}: pack {index}'s '{limit}' claims give no strength at \
                {:.4} mm (missing data): {e}", mm(d)),
            None => format!("string material pack {index}: {e}"),
        })?;
    let property = state.property(limit).ok_or("tensile limit resolution lost its property")?;
    let receipt = &property.answer().receipt;
    // The selected claim's own rows: the bracket a linear rule interpolated, or the row it hit.
    let curve = match pack.card().claims().claim(receipt.selected).map(|claim| &claim.value) {
        Some(PropertyValue::Curve { abscissa, knots, .. }) => Some((abscissa.as_str(), knots.as_slice())),
        _ => None,
    };
    let row = |x: f64| curve.and_then(|(_, knots)| knots.iter().find(|knot| knot.0 == x).copied())
        .ok_or_else(|| format!("string material key {key}: the '{limit}' answer lost its tabulated row at {x}"));
    let axis = curve.map_or("", |(abscissa, _)| abscissa);
    let (decision, rows) = match &receipt.decision {
        EvaluationDecision::ExactTabulated { at } => (format!("tabulated {axis} row at {at}"), vec![row(*at)?]),
        EvaluationDecision::LinearInside { x_lo, x_hi } =>
            (format!("linear in {axis} between rows at {x_lo} and {x_hi}"), vec![row(*x_lo)?, row(*x_hi)?]),
        EvaluationDecision::ConstantWithinValidity => ("constant over its validity".to_string(), vec![]),
        EvaluationDecision::ExactScalar => ("single tabulated value".to_string(), vec![]),
        other => (format!("{other:?}"), vec![]),
    };
    let caveat = match (domain, query_at) {
        (None, _) => Some(format!("strength not diameter-resolved; '{limit}' from {} states no wire diameter, \
            diameter effect unbounded", pack.pack_id())),
        (Some((lo, hi)), Some(at)) if at != diameter_m => Some(if lo == hi {
            format!("strength extrapolated; '{limit}' from {} holds only for {:.4} mm wire and is taken at every \
                core diameter, diameter effect unbounded", pack.pack_id(), mm(lo))
        } else {
            format!("strength extrapolated; '{limit}' from {} holds for {:.4}-{:.4} mm wire; a core outside takes \
                the nearest stated diameter's value", pack.pack_id(), mm(lo), mm(hi))
        }),
        _ => None,
    };
    Ok((property.clone(), LimitReceipt { pack: index, diameter_m, evaluated_at_m: query_at, domain_m: domain,
        limit_pa: property.value_si(), decision, rows, extrapolated: query_at != Some(diameter_m), caveat }))
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
    /// Per key, where each core layer's limit came from (spatial layer order).
    limit_receipts: BTreeMap<u8, Vec<LimitReceipt>>,
    tension: TensionPolicy,
    limit_property: String,
    non_physical_override: bool,
    extrapolate_diameter: bool,
    diameter_free: bool,
    strength_pack: Option<usize>,
    /// Keys whose core and wrap diameters the construction declares back-solved.
    back_solved: Option<(u8, u8)>,
    packs: Vec<Arc<NormalizedMaterialCardPack>>,
}

#[derive(Clone, Copy)]
enum Role { Core, Wrap, Strength }

struct Declaration {
    tension: TensionPolicy,
    loss: BendingLossLaw,
    limit_property: String,
    non_physical_override: bool,
    extrapolate_diameter: bool,
    diameter_free: bool,
    strength_pack: Option<usize>,
    back_solved: Option<(u8, u8)>,
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
    let strength = lines.next().map(|line| line.split(',').map(str::trim).collect::<Vec<_>>());
    let Some(["strength", name, flags @ ..]) = strength.as_deref() else {
        return Err("third line must be strength,<tensile-limit property>[,non-physical-override]\
            [,extrapolate-diameter][,diameter-free][,pack=<index>]".into());
    };
    let limit_property = property_name(name, "tensile-limit")?;
    let (mut non_physical_override, mut extrapolate_diameter, mut strength_pack) = (false, false, None);
    let mut diameter_free = false;
    for &flag in flags {
        match flag {
            "non-physical-override" if !non_physical_override => non_physical_override = true,
            "extrapolate-diameter" if !extrapolate_diameter => extrapolate_diameter = true,
            "diameter-free" if !diameter_free => diameter_free = true,
            _ if flag.starts_with("pack=") && strength_pack.is_none() => {
                let index = flag["pack=".len()..].parse::<usize>().map_err(|_| "invalid strength pack index")?;
                if index >= packs { return Err("strength pack index is out of range".into()); }
                strength_pack = Some(index);
            }
            _ => return Err(format!("unknown or repeated strength flag '{}': flags are non-physical-override, \
                extrapolate-diameter, diameter-free and pack=<index>", escape(flag))),
        }
    }
    let mut lines = lines.peekable();
    let loss = match lines.peek().map(|line| line.split(',').map(str::trim).collect::<Vec<_>>()).as_deref() {
        Some(["loss", "kelvin-voigt"]) => { lines.next(); BendingLossLaw::KelvinVoigt }
        Some(["loss", "loss-angle"]) => { lines.next(); BendingLossLaw::LossAngle }
        Some(["loss", "thermoelastic"]) => { lines.next(); BendingLossLaw::Thermoelastic }
        Some(["loss", ..]) => return Err("loss line must be loss,kelvin-voigt, loss,loss-angle or loss,thermoelastic".into()),
        _ => BendingLossLaw::KelvinVoigt,
    };
    let back_solved = match lines.peek().map(|line| line.split(',').map(str::trim).collect::<Vec<_>>()).as_deref() {
        Some(["geometry", "back-solved", keys]) => {
            let range = keys.split_once('-').and_then(|(lo, hi)| Some((lo.parse::<u8>().ok()?, hi.parse::<u8>().ok()?)));
            let Some((lo, hi)) = range.filter(|&(lo, hi)| (21..=108).contains(&lo) && (lo..=108).contains(&hi)) else {
                return Err("geometry,back-solved needs a key range <first>-<last> within 21-108".into());
            };
            lines.next();
            Some((lo, hi))
        }
        Some(["geometry", ..]) => return Err("geometry line must be geometry,back-solved,<first>-<last>".into()),
        _ => None,
    };
    if lines.next() != Some(COLUMNS) {
        return Err(format!("the line after strength (and loss, geometry) must be {COLUMNS}"));
    }
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
    if let Some((lo, hi)) = back_solved
        && let Some(key) = (lo..=hi).find(|key| !cores.contains_key(key)) {
        return Err(format!("geometry,back-solved names key {key}, which has no section"));
    }
    Ok(Declaration { tension, loss, limit_property, non_physical_override, extrapolate_diameter, diameter_free,
        strength_pack, back_solved,
        cores, wraps })
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
        if let Some(index) = declaration.strength_pack { roles[index].push(Role::Strength); }
        let mut packs = Vec::with_capacity(bytes.len());
        let mut states = Vec::with_capacity(bytes.len());
        let mut loss_rated = Vec::with_capacity(bytes.len());
        let mut excess_measured = Vec::with_capacity(bytes.len());
        for (index, (bytes, roles)) in bytes.iter().zip(&roles).enumerate() {
            let pack = Arc::new(NormalizedMaterialCardPack::from_bytes(bytes)
                .map_err(|e| e.to_string())?);
            let core = roles.iter().any(|role| matches!(role, Role::Core));
            if roles.is_empty() { return Err(format!("string material pack {index} is not used by any layer")); }
            // A pack used only for strength holds no section state: its limit is
            // evaluated per core diameter in `from_states`.
            if roles.iter().all(|role| matches!(role, Role::Strength)) {
                packs.push(pack);
                states.push(None);
                loss_rated.push(true);
                excess_measured.push(true);
                continue;
            }
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
            let refuse = |e: MaterialStatePointError| format!("string material pack {index}: {e}");
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
            states.push(Some(state));
            loss_rated.push(rated);
            excess_measured.push(!dropped.iter().any(|property| property == viscosity_property));
        }
        Self::from_states(declaration, &states, &loss_rated, &excess_measured, viscosity_property, packs, point)
    }

    #[allow(clippy::too_many_arguments)]
    fn from_states(declaration: Declaration, states: &[Option<ResolvedMaterialStatePoint>], loss_rated: &[bool],
        excess_measured: &[bool], viscosity_property: &str, packs: Vec<Arc<NormalizedMaterialCardPack>>,
        point: &QueryPoint) -> Result<Self, String> {
        let Declaration { tension, loss, limit_property, non_physical_override, extrapolate_diameter, diameter_free,
            strength_pack, back_solved, cores, mut wraps } = declaration;
        let flags = DiameterFlags { extrapolate: extrapolate_diameter, diameter_free };
        let state = |index: usize| states[index].clone()
            .ok_or_else(|| format!("string material pack {index} is a strength source only and holds no section"));
        let mut limit_receipts = BTreeMap::new();
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
            let sources = layers.iter().map(|&(_, index)| strength_pack.unwrap_or(index)).collect::<Vec<_>>();
            let layers = layers.into_iter().map(|(outer_radius_m, index)|
                Ok(ConcentricMaterialLayer { outer_radius_m, material: state(index)? }))
                .collect::<Result<Vec<_>, String>>()?;
            let elastic = ConcentricElasticSection::try_new(0.0, layers).map_err(failed)?;
            // The drawn wire's diameter: every layer's claim is evaluated there.
            let diameter_m = 2.0 * elastic.outer_radius_m();
            let (limits, receipts): (Vec<ResolvedScalarProperty>, Vec<LimitReceipt>) = sources.iter()
                .map(|&source| resolve_limit(&packs[source], point, &limit_property, source, key, diameter_m, flags))
                .collect::<Result<Vec<_>, String>>()?.into_iter().unzip();
            // A limit that states no diameter bounds a check but cannot set a playing tension.
            if let (TensionPolicy::YieldFraction(_), Some(free)) =
                (tension, receipts.iter().find(|receipt| receipt.domain_m.is_none())) {
                return Err(format!("string material key {key}: tension,yield-fraction sets the playing tension \
                    from the limit, but '{limit_property}' from {} states no wire diameter (strength not \
                    diameter-resolved) and no sourced yield stress is declared; nothing back-solves a yield \
                    (missing data)", packs[free.pack].pack_id()));
            }
            capacities.insert(key, elastic.axial_tensile_capacity_from_resolved(&limit_property, &limits)
                .map_err(failed)?);
            limit_receipts.insert(key, receipts);
            if let Some(rows) = wraps.remove(&key) {
                let mut under = 2.0 * elastic.outer_radius_m();
                let mut helical = Vec::with_capacity(rows.len());
                for (wire_diameter_m, pitch, index) in rows {
                    let pitch_m = match pitch {
                        Some(pitch) => pitch,
                        None => close_wound_pitch_m(under, wire_diameter_m).map_err(failed)?,
                    };
                    helical.push(HelicalWrap { wire_diameter_m, pitch_m, material: state(index)? });
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
            capacities, limit_receipts, tension, limit_property, non_physical_override, extrapolate_diameter,
            diameter_free, strength_pack, back_solved, packs })
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
    /// The construction declared `extrapolate-diameter`.
    pub fn extrapolate_diameter(&self) -> bool { self.extrapolate_diameter }
    /// The construction declared `diameter-free`.
    pub fn diameter_free(&self) -> bool { self.diameter_free }
    /// Keys whose core and wrap diameters the construction declares back-solved.
    pub fn back_solved_keys(&self) -> Option<(u8, u8)> { self.back_solved }
    /// What a declared back-solved geometry means for every result derived from it.
    pub fn geometry_caveat(&self) -> Option<String> {
        self.back_solved.map(|(lo, hi)| format!("core and wrap diameters of keys {lo}-{hi} are back-solved, not \
            measured or cited; every tension, stress, limit, pitch and loss computed from them inherits that"))
    }
    /// The pack named by `pack=<index>` as every core's strength source.
    pub fn strength_pack(&self) -> Option<usize> { self.strength_pack }
    /// Per key, where each core layer's tensile limit came from.
    pub fn limit_receipts(&self) -> &BTreeMap<u8, Vec<LimitReceipt>> { &self.limit_receipts }
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
            let receipt = &self.limit_receipts[&course.midi][layer];
            keys.push(KeyStress {
                key: course.midi, governing_layer: layer, max_member_tension_n: tension,
                required_stress_pa: capacity.layer_stress_pa(layer, tension).map_err(fail)?,
                limit_pa: capacity.layer_limits_pa()[layer],
                utilization: capacity.utilization(tension).map_err(fail)?,
                capacity_n: capacity.capacity_n(), f_max_hz,
                limit_evaluated_at_m: receipt.evaluated_at_m, limit_extrapolated: receipt.extrapolated,
            });
        }
        let mut limit_caveats: Vec<String> = Vec::new();
        for caveat in courses.iter().flat_map(|course| &self.limit_receipts[&course.midi])
            .filter_map(|receipt| receipt.caveat.as_ref()) {
            if !limit_caveats.contains(caveat) { limit_caveats.push(caveat.clone()); }
        }
        Ok(StrengthReport { limit_property: self.limit_property.clone(),
            non_physical_override: self.non_physical_override, keys, limit_caveats,
            geometry_caveat: self.geometry_caveat() })
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
            let limits = self.limit_receipts[key].iter().map(LimitReceipt::json).collect::<Vec<_>>().join(",");
            format!(concat!("{{\"key\":{},\"stateIdentities\":[{}],\"outerRadiiM\":{:?},",
                "\"lineDensityKgM\":{},\"axialRigidityN\":{},\"flexuralRigidityNm2\":{},",
                "\"viscousRigidityNm2S\":{},\"hystereticRigidityNm2\":{},\"omegaBandRadS\":{},\"thermoelastic\":{},\"wraps\":{},",
                "\"outerRadiusM\":{},\"tensileCapacityN\":{},\"governingLayer\":{},\"layerLimitsPa\":{:?},",
                "\"limitReceipts\":[{}]}}"),
                key, identities, radii, elastic.line_density_kg_m()
                    + self.wrappings.get(key).map_or(0.0, HelicalWrapping::line_density_kg_m),
                elastic.axial_stiffness_n(), elastic.bending_stiffness_n_m2(),
                viscous, hysteretic, band, thermoelastic, wraps,
                self.outer_radius_m(*key).unwrap_or(0.0), capacity.capacity_n(),
                capacity.governing_layer(), capacity.layer_limits_pa(), limits)
        }).collect::<Vec<_>>().join(",");
        let receipt = match self.loss {
            BendingLossLaw::Thermoelastic => "thermoelastic loss is computed at every mode; a measured anelastic \
                excess is added only inside its source band; elsewhere intrinsic anelastic loss is omitted, so the \
                decay there is a lower bound",
            _ => "measured bending loss only; a mode outside its source band refuses",
        };
        format!(concat!("{{\"model\":\"bonded-core-kelvin-voigt-helical-wraps\",\"importedPacks\":{},",
            "\"tensionPolicy\":{},\"bendingLoss\":\"{}\",\"lossReceipt\":\"{}\",\"missingLossData\":{},",
            "\"limitProperty\":\"{}\",\"nonPhysicalOverride\":{},\"extrapolateDiameter\":{},\"diameterFree\":{},",
            "\"strengthPack\":{},\"geometry\":{},",
            "\"sections\":[{}]}}"),
            self.packs.len(), self.tension.json(), self.loss.name(), receipt, self.missing_loss.is_some(),
            escape(&self.limit_property), self.non_physical_override, self.extrapolate_diameter, self.diameter_free,
            self.strength_pack.map_or_else(|| "null".to_string(), |index| index.to_string()),
            self.back_solved.zip(self.geometry_caveat()).map_or_else(|| "null".to_string(), |((lo, hi), caveat)|
                format!("{{\"backSolvedKeys\":[{lo},{hi}],\"caveat\":\"{}\"}}", escape(&caveat))), sections)
    }
}

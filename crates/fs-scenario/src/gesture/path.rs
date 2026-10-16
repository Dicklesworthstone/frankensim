//! Typed writable paths and the isolated v1 alias vocabulary.
//!
//! Writable declarations are the union of the description's configuration
//! assignments. An absent parameter is never inferred from a gesture. Legacy
//! indices select matching declarations in part-declaration order (root first),
//! not lexical display-name order. No geometry or mechanism is synthesized.
use super::{GestureError, GestureSchedule, GestureTarget, GestureValue};
use crate::assembly::{
    AssemblyDescription, DEFAULT_ASSEMBLY_BUDGET, ParameterKey, ParameterPath, PartKind,
};
use crate::entity::{EntityId, EntityKind};
use fs_blake3::ContentHash;
use fs_qty::Dims;
use std::collections::BTreeMap;

pub const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
pub const FORCE: Dims = Dims([1, 1, -2, 0, 0, 0]);
pub const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);
pub const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);

pub fn path_name(path: &ParameterPath) -> String {
    format!("{}:{}", path.entity.token(), path.parameter.name)
}
pub fn value_dimensions(value: &GestureValue) -> Dims {
    match value {
        GestureValue::LengthM(_) => LENGTH,
        GestureValue::TensionN(_) => FORCE,
        GestureValue::PressurePa(_) => PRESSURE,
        GestureValue::VelocityMPerS(_) | GestureValue::StrikeVelocity { .. } => SPEED,
        // Composite commands retain their component units in GestureValue.
        GestureValue::Fraction(_)
        | GestureValue::Multiplier(_)
        | GestureValue::AngleRad(_)
        | GestureValue::Fret { .. }
        | GestureValue::Bow { .. }
        | GestureValue::TerminationSwap { .. } => Dims::NONE,
    }
}
pub(crate) fn accepts(path: &ParameterPath, value: &GestureValue) -> bool {
    if path.parameter.dims != value_dimensions(value) {
        return false;
    }
    match path.parameter.name.as_str() {
        "inserted_length" | "slide_length" | "rest_aperture" => {
            matches!(value, GestureValue::LengthM(_))
        }
        "open_fraction" | "sustain" | "una_corda" | "sostenuto_hold" => {
            matches!(value, GestureValue::Fraction(_))
        }
        "tension" => matches!(value, GestureValue::TensionN(_)),
        "blowing_pressure" => matches!(value, GestureValue::PressurePa(_)),
        "prestress" => matches!(value, GestureValue::Multiplier(_)),
        "jet_speed" => matches!(value, GestureValue::VelocityMPerS(_)),
        "jet_angle" => matches!(value, GestureValue::AngleRad(_)),
        "fret" => matches!(value, GestureValue::Fret { .. }),
        "bow" => matches!(value, GestureValue::Bow { .. }),
        "strike" => matches!(value, GestureValue::StrikeVelocity { .. }),
        "termination" => matches!(value, GestureValue::TerminationSwap { .. }),
        _ => true,
    }
}
pub(crate) fn validate_path(path: &ParameterPath) -> Result<(), GestureError> {
    if path.parameter.name.is_empty()
        || path.parameter.name.chars().any(char::is_whitespace)
        || path.parameter.name.chars().any(char::is_control)
    {
        return Err(GestureError::UnknownParameter {
            path: path_name(path),
        });
    }
    Ok(())
}
pub(crate) fn encode(path: &ParameterPath) -> String {
    let dims = path
        .parameter
        .dims
        .0
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join("\t");
    format!(
        "param\t{}\t{}\t{dims}",
        path.entity.token(),
        path.parameter.name
    )
}
pub(crate) fn decode(s: &str) -> Option<ParameterPath> {
    let mut c = s.split('\t');
    if c.next()? != "param" {
        return None;
    }
    let (kind, hex) = c.next()?.split_once(':')?;
    let kind = match kind {
        "assembly" => EntityKind::Assembly,
        "part" => EntityKind::Part,
        "region" => EntityKind::Region,
        "surface" => EntityKind::Surface,
        "interface" => EntityKind::Interface,
        _ => return None,
    };
    if hex.len() != 64
        || hex
            .bytes()
            .any(|v| !v.is_ascii_digit() && !(b'a'..=b'f').contains(&v))
    {
        return None;
    }
    let entity = EntityId::from_wire(kind, ContentHash::from_hex(hex)?);
    let name = c.next()?.to_string();
    let mut dims = [0i8; 6];
    for d in &mut dims {
        *d = c.next()?.parse().ok()?;
    }
    if c.next().is_some() {
        return None;
    }
    let path = ParameterPath {
        entity,
        parameter: ParameterKey::new(name, Dims(dims)),
    };
    validate_path(&path).ok()?;
    Some(path)
}

/// The only module interpreting the legacy u8 member indices. New consumers
/// bind ParameterPath; this descriptor exists only for compatibility adapters.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LegacyAlias {
    pub parameter: ParameterKey,
    pub index: usize,
    pub part_kind: Option<PartKind>,
}
pub fn legacy_alias(target: &GestureTarget) -> Option<LegacyAlias> {
    let (name, dims, index) = match target {
        GestureTarget::ValveInsertedLength { valve } => {
            ("inserted_length", LENGTH, usize::from(*valve))
        }
        GestureTarget::ToneHoleSigma { hole } => ("open_fraction", Dims::NONE, usize::from(*hole)),
        GestureTarget::SlideLength => ("slide_length", LENGTH, 0),
        GestureTarget::StringTension { string } => ("tension", FORCE, usize::from(*string)),
        GestureTarget::Fret { string } => ("fret", Dims::NONE, usize::from(*string)),
        GestureTarget::BowStroke { string } => ("bow", Dims::NONE, usize::from(*string)),
        GestureTarget::HammerStrike { string } => ("strike", SPEED, usize::from(*string)),
        GestureTarget::SustainPedal => ("sustain", Dims::NONE, 0),
        GestureTarget::UnaCordaPedal => ("una_corda", Dims::NONE, 0),
        GestureTarget::SostenutoHold { string } => {
            ("sostenuto_hold", Dims::NONE, usize::from(*string))
        }
        GestureTarget::TerminationSwap => ("termination", Dims::NONE, 0),
        GestureTarget::BlowingPressure => ("blowing_pressure", PRESSURE, 0),
        GestureTarget::PreStress => ("prestress", Dims::NONE, 0),
        GestureTarget::RestAperture => ("rest_aperture", LENGTH, 0),
        GestureTarget::JetSpeed => ("jet_speed", SPEED, 0),
        GestureTarget::JetAngle => ("jet_angle", Dims::NONE, 0),
        GestureTarget::Parameter(_) => return None,
    };
    Some(LegacyAlias {
        parameter: ParameterKey::new(name, dims),
        index,
        part_kind: match name {
            "inserted_length" => Some(PartKind::ValvePattern),
            "open_fraction" => Some(PartKind::Duct),
            "tension" | "fret" | "bow" | "strike" | "sostenuto_hold" => Some(PartKind::String),
            _ => None,
        },
    })
}
/// Decode compatibility selection without exposing u8 matching to physics.
pub fn legacy_member(target: &GestureTarget, key: &str) -> Option<usize> {
    let a = legacy_alias(target)?;
    (a.parameter.name == key).then_some(a.index)
}
/// A schedule whose every path was checked against one immutable description.
#[derive(Clone, Debug)]
pub struct ResolvedSchedule {
    description: ContentHash,
    source: ContentHash,
    schedule: GestureSchedule,
}
impl ResolvedSchedule {
    pub fn description_hash(&self) -> ContentHash {
        self.description
    }
    /// Original v1/v2 identity for the explicit alias-migration receipt.
    pub fn source_hash(&self) -> ContentHash {
        self.source
    }
    pub fn schedule(&self) -> &GestureSchedule {
        &self.schedule
    }
    pub fn track(&self, path: &ParameterPath) -> Result<&super::GestureTrack, GestureError> {
        self.schedule
            .tracks()
            .iter()
            .find(|t| matches!(&t.target,GestureTarget::Parameter(p) if p==path))
            .ok_or_else(|| GestureError::UnknownParameter {
                path: path_name(path),
            })
    }
}
pub fn declared_parameters(
    description: &AssemblyDescription,
) -> Result<Vec<ParameterPath>, GestureError> {
    let mut parameters = BTreeMap::<(EntityId, String), Dims>::new();
    for state in &description.input().configurations {
        for (p, _) in &state.assignments {
            validate_path(p)?;
            if let Some(dims) =
                parameters.insert((p.entity, p.parameter.name.clone()), p.parameter.dims)
            {
                if dims != p.parameter.dims {
                    return Err(GestureError::DimensionMismatch {
                        path: path_name(p),
                        expected: dims,
                        actual: p.parameter.dims,
                    });
                }
            }
        }
    }
    let mut paths: Vec<_> = parameters
        .into_iter()
        .map(|((entity, name), dims)| ParameterPath {
            entity,
            parameter: ParameterKey::new(name, dims),
        })
        .collect();
    let rank = |id| {
        if id == description.input().root {
            0
        } else {
            description
                .input()
                .parts
                .iter()
                .position(|p| p.entity == id)
                .map_or(usize::MAX, |i| i + 1)
        }
    };
    paths.sort_by_key(|p| (rank(p.entity), p.entity.token(), p.parameter.name.clone()));
    Ok(paths)
}
pub fn resolve_aliases(
    schedule: &GestureSchedule,
    description: &AssemblyDescription,
) -> Result<ResolvedSchedule, GestureError> {
    let schedule_source_hash = schedule.content_hash();
    let declarations = declared_parameters(description)?;
    let mut tracks = schedule.tracks().to_vec();
    for t in &mut tracks {
        let path = match &t.target {
            GestureTarget::Parameter(p) => {
                let declared = declarations
                    .iter()
                    .find(|d| d.entity == p.entity && d.parameter.name == p.parameter.name)
                    .ok_or_else(|| GestureError::UnknownParameter { path: path_name(p) })?;
                if declared.parameter.dims != p.parameter.dims {
                    return Err(GestureError::DimensionMismatch {
                        path: path_name(p),
                        expected: declared.parameter.dims,
                        actual: p.parameter.dims,
                    });
                }
                p.clone()
            }
            alias => {
                let a = legacy_alias(alias).ok_or(GestureError::Invalid {
                    what: "unrecognized alias",
                })?;
                let resolved = declarations
                    .iter()
                    .filter(|p| {
                        p.parameter == a.parameter
                            && a.part_kind.is_none_or(|kind| {
                                description
                                    .input()
                                    .parts
                                    .iter()
                                    .any(|part| part.entity == p.entity && part.kind == kind)
                            })
                    })
                    .nth(a.index)
                    .cloned();
                resolved.ok_or(GestureError::AliasUnresolved {
                    target: a.parameter.name,
                    index: a.index,
                })?
            }
        };
        t.target = GestureTarget::Parameter(path);
    }
    let schedule = GestureSchedule::try_new_v2(schedule.control_rate_hz, tracks)?;
    let description = description
        .content_hash(DEFAULT_ASSEMBLY_BUDGET)
        .map_err(|_| GestureError::Invalid {
            what: "description hash budget",
        })?;
    Ok(ResolvedSchedule {
        description,
        source: schedule_source_hash,
        schedule,
    })
}

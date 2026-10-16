//! Pre-render admission, using the existing registry and physics owners.
//! Numerical domain policies are caller-declared evidence bounds, not new
//! registry rows or a proof of specimen fidelity. No band is silently narrowed.
use crate::material_assembly::MaterialSource;
use crate::selector::{Gate, RegistryRow, RegistrySnapshot};
use fs_blake3::ContentHash;
use fs_duct::{LossModel, Segment};
use fs_evidence::ValidityDomain;
use fs_material::gas::GasState;
use fs_material::state_point::{
    DENSITY_PROPERTY, ScalarAdmissibility, ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY,
    resolve_material_state_point,
};
use fs_qty::{Dims, QuantitySpec};
use fs_scenario::assembly::{AssemblyDescription, GeometrySource, PartDecl, PartKind};
use fs_scenario::{BeatingReed, EntityId};
use std::collections::BTreeSet;

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const HZ: Dims = Dims([0, 0, -1, 0, 0, 0]);
const TAU: f64 = core::f64::consts::TAU;

/// Conservative command/response envelope, all SI. The lowest excited frequency
/// is explicit even when the requested observation band includes DC. Response
/// displacement must come from the declared experiment; force alone cannot
/// determine a string's Keulegan–Carpenter number.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DynamicsEnvelope {
    pub max_mouth_pressure_pa: f64,
    pub max_bow_force_n: f64,
    pub max_pluck_force_n: f64,
    pub lowest_excitation_hz: f64,
    pub max_string_displacement_m: f64,
}

/// Exact identity in the one existing instrument-claims registry.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RegistryClaim {
    pub filling: String,
    pub image: String,
    pub qoi: String,
}

/// Selected generic physics owner. This selection never follows a part's name.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DescriptionModel {
    CylindricalDuct { loss: LossModel },
    CompactToneHole,
    CompactRadiation { flange: fs_phs::MouthFlange },
    LumpedReed,
    OscillatingCylinderString,
}
impl DescriptionModel {
    pub const fn module(self) -> &'static str {
        match self {
            Self::CylindricalDuct { .. } => "fs-duct/cylindrical",
            Self::CompactToneHole => "fs-duct/compact-tone-hole",
            Self::CompactRadiation { .. } => "fs-phs/compact-radiation",
            Self::LumpedReed => "fs-couple/lumped-reed",
            Self::OscillatingCylinderString => "fs-couple/oscillating-cylinder-string",
        }
    }
}

/// One selected part/model and its existing typed validity box. Multiple owners
/// can act on one part, but duplicate (part,module) entries refuse. Required
/// policy axes are never unconstrained: frequency_hz for all models, ka for
/// radiation, edge_displacement_ratio for holes, and length_over_shock_distance
/// plus evanescent_relative_error for ducts. Extra axes must be supplied in the
/// description's dimensioned parameter table as well.
#[derive(Clone, Debug)]
pub struct ModuleDomain<'a> {
    pub part: EntityId,
    pub model: DescriptionModel,
    pub claim: RegistryClaim,
    pub domain: ValidityDomain,
    /// Exact resolved card/query bindings, one per description material ref.
    pub materials: Vec<MaterialSource<'a>>,
}

/// All inputs to admission are explicit. RegistrySnapshot is the same type the
/// runtime selector consumes; this module has no parallel claims catalogue.
#[derive(Clone, Debug)]
pub struct ClaimDomain<'a> {
    pub band_hz: (f64, f64),
    pub dynamics: DynamicsEnvelope,
    pub gas: GasState,
    pub registry: &'a RegistrySnapshot,
    pub modules: Vec<ModuleDomain<'a>>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct DomainAxis {
    pub module: &'static str,
    pub axis: String,
    pub value: f64,
    pub domain: (f64, f64),
    /// Owner geometry constraints can have a strict upper endpoint.
    pub upper_inclusive: bool,
    pub unit: Dims,
    pub refinement: Option<&'static str>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct ModelFormLimitation {
    pub part: EntityId,
    pub module: &'static str,
    /// None is the construction state; every named configuration is checked too.
    pub configuration: Option<String>,
    pub statement: &'static str,
    pub claim: RegistryRow,
    pub axes: Vec<DomainAxis>,
    pub material_states: Vec<ContentHash>,
}

/// Only validate can construct this token. It borrows the exact immutable
/// description and claim: description render consumers must accept this token,
/// rather than a bypass bool or an unrelated second description. Admission
/// checks model domains, not Q3/Q4
/// qualification; even a green registry row does not establish specimen fidelity.
#[derive(Debug)]
pub struct Admission<'a> {
    description: &'a AssemblyDescription,
    claim: &'a ClaimDomain<'a>,
    limitations: Vec<ModelFormLimitation>,
}
impl<'a> Admission<'a> {
    pub fn description(&self) -> &'a AssemblyDescription {
        self.description
    }
    pub fn claim(&self) -> &'a ClaimDomain<'a> {
        self.claim
    }
    pub fn claim_band_hz(&self) -> (f64, f64) {
        self.claim.band_hz
    }
    pub fn limitations(&self) -> &[ModelFormLimitation] {
        &self.limitations
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum ValidityError {
    OutOfDomain {
        part: EntityId,
        axis: DomainAxis,
    },
    MissingAxis {
        part: EntityId,
        module: &'static str,
        axis: String,
    },
    DimensionMismatch {
        part: EntityId,
        module: &'static str,
        axis: String,
        expected: Dims,
        found: Dims,
    },
    UnknownRequiredProperty {
        part: EntityId,
        module: &'static str,
        property: String,
    },
    UnsupportedPart {
        part: EntityId,
        kind: PartKind,
    },
    InvalidBinding {
        part: EntityId,
        module: &'static str,
        reason: String,
    },
    RegistryRefusal {
        part: EntityId,
        claim: RegistryClaim,
    },
    OwnerRefusal {
        part: EntityId,
        module: &'static str,
        reason: String,
    },
}
impl std::fmt::Display for ValidityError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "assembly validity: {self:?}")
    }
}
impl std::error::Error for ValidityError {}

fn missing(part: &PartDecl, module: &'static str, axis: &str) -> ValidityError {
    ValidityError::MissingAxis {
        part: part.entity,
        module,
        axis: axis.into(),
    }
}
fn binding(part: EntityId, module: &'static str, reason: impl Into<String>) -> ValidityError {
    ValidityError::InvalidBinding {
        part,
        module,
        reason: reason.into(),
    }
}
fn owner(part: &PartDecl, module: &'static str, error: impl std::fmt::Display) -> ValidityError {
    ValidityError::OwnerRefusal {
        part: part.entity,
        module,
        reason: error.to_string(),
    }
}
fn params(
    part: &PartDecl,
) -> Result<
    (
        &str,
        &[(fs_scenario::assembly::ParameterKey, fs_qty::QtyAny)],
    ),
    ValidityError,
> {
    match &part.geometry {
        GeometrySource::Parametric { generator, params } => Ok((generator, params)),
        _ => Err(ValidityError::UnsupportedPart {
            part: part.entity,
            kind: part.kind,
        }),
    }
}
fn parameter(
    part: &PartDecl,
    module: &'static str,
    name: &str,
    dims: Dims,
) -> Result<f64, ValidityError> {
    let (_, table) = params(part)?;
    let (key, value) = table
        .iter()
        .find(|(key, _)| key.name == name)
        .ok_or_else(|| missing(part, module, name))?;
    if key.dims != dims || value.dims != dims {
        return Err(ValidityError::DimensionMismatch {
            part: part.entity,
            module,
            axis: name.into(),
            expected: dims,
            found: value.dims,
        });
    }
    Ok(value.value)
}
fn geometry(part: &PartDecl, module: &'static str, expected: &str) -> Result<(), ValidityError> {
    let (name, _) = params(part)?;
    // A fixture marker changes provenance, never geometry semantics.
    let name = match name.strip_prefix("test-fixture:") {
        Some(name) => name,
        None => name,
    };
    if name != expected {
        return Err(ValidityError::OutOfDomain {
            part: part.entity,
            axis: DomainAxis {
                module,
                axis: "circular_cross_section".into(),
                value: 0.0,
                domain: (1.0, 1.0),
                upper_inclusive: true,
                unit: Dims::NONE,
                refinement: Some("M2.4-W1c / M3.2"),
            },
        });
    }
    Ok(())
}
fn axis(
    part: &PartDecl,
    module: &'static str,
    name: &str,
    value: f64,
    domain: (f64, f64),
    unit: Dims,
    refinement: Option<&'static str>,
) -> Result<DomainAxis, ValidityError> {
    let a = DomainAxis {
        module,
        axis: name.into(),
        value,
        domain,
        upper_inclusive: true,
        unit,
        refinement,
    };
    if !value.is_finite()
        || !domain.0.is_finite()
        || !domain.1.is_finite()
        || domain.0 > domain.1
        || value < domain.0
        || value > domain.1
    {
        Err(ValidityError::OutOfDomain {
            part: part.entity,
            axis: a,
        })
    } else {
        Ok(a)
    }
}
fn policy(
    part: &PartDecl,
    m: &ModuleDomain<'_>,
    name: &str,
    value: f64,
    unit: Dims,
    refinement: Option<&'static str>,
) -> Result<DomainAxis, ValidityError> {
    let domain = m
        .domain
        .bound(name)
        .ok_or_else(|| missing(part, m.model.module(), name))?;
    if let Some(quantity) = m.domain.axis_quantities().get(name)
        && *quantity != QuantitySpec::dimensional(unit)
    {
        return Err(ValidityError::DimensionMismatch {
            part: part.entity,
            module: m.model.module(),
            axis: name.into(),
            expected: unit,
            found: quantity.dims(),
        });
    }
    axis(
        part,
        m.model.module(),
        name,
        value,
        domain,
        unit,
        refinement,
    )
}

fn materials(part: &PartDecl, m: &ModuleDomain<'_>) -> Result<Vec<ContentHash>, ValidityError> {
    let module = m.model.module();
    let requirements = if matches!(
        m.model,
        DescriptionModel::LumpedReed | DescriptionModel::OscillatingCylinderString
    ) {
        vec![
            ScalarPropertyRequirement::try_new(
                YOUNG_MODULUS_PROPERTY,
                PRESSURE,
                ScalarAdmissibility::StrictlyPositive,
            ),
            ScalarPropertyRequirement::try_new(
                DENSITY_PROPERTY,
                Dims([-3, 1, 0, 0, 0, 0]),
                ScalarAdmissibility::StrictlyPositive,
            ),
        ]
        .into_iter()
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| owner(part, module, e))?
    } else {
        vec![]
    };
    // Explicit Unknown has priority over absent card data, including when the
    // author has not yet provided the resolver binding.
    for material in &part.materials {
        for r in &requirements {
            if material.require_known(r.name()).is_err() {
                return Err(ValidityError::UnknownRequiredProperty {
                    part: part.entity,
                    module,
                    property: r.name().into(),
                });
            }
        }
    }
    if (!requirements.is_empty() && part.materials.is_empty())
        || part.materials.len() != m.materials.len()
    {
        return Err(binding(
            part.entity,
            module,
            "exact material sources required for description refs",
        ));
    }
    let mut states = Vec::new();
    for (reference, source) in part.materials.iter().zip(&m.materials) {
        if reference.card != source.card.content_hash()
            || reference.query_point.len() != source.point.axes().len()
        {
            return Err(binding(
                part.entity,
                module,
                "material card/query does not match description",
            ));
        }
        for (key, value) in &reference.query_point {
            if source.point.axes().get(&key.name) != Some(&value.value)
                || source.point.axis_quantities().get(&key.name)
                    != Some(&QuantitySpec::dimensional(key.dims))
            {
                return Err(binding(
                    part.entity,
                    module,
                    "material query quantity/value does not match description",
                ));
            }
        }
        if !requirements.is_empty() {
            let state = resolve_material_state_point(
                source.card,
                source.point,
                &requirements,
                source.selection.clone(),
            )
            .map_err(|e| owner(part, module, e))?;
            states.push(state.identity());
        }
    }
    Ok(states)
}

/// G0/G1 model-domain admission, before constructing any render state.
/// Unsupported geometry, missing axes/materials, nonfinite values, unknown or
/// refused registry rows and all exceeded bounds fail closed. No render or fit
/// is performed and no source coefficient is tuned here.
pub fn validate<'a>(
    description: &'a AssemblyDescription,
    claim: &'a ClaimDomain<'a>,
) -> Result<Admission<'a>, ValidityError> {
    let input = description.input();
    let root = input.root;
    crate::acoustic_realize::validate_external_gas(&claim.gas)
        .map_err(|e| binding(root, "fs-material/gas", e.to_string()))?;
    let (lo, hi) = claim.band_hz;
    if !lo.is_finite() || !hi.is_finite() || lo < 0.0 || hi <= 0.0 || lo > hi {
        return Err(binding(
            root,
            "claim",
            "finite ordered nonnegative frequency band required",
        ));
    }
    let d = claim.dynamics;
    if ![
        d.max_mouth_pressure_pa,
        d.max_bow_force_n,
        d.max_pluck_force_n,
        d.max_string_displacement_m,
    ]
    .into_iter()
    .all(|v| v.is_finite() && v >= 0.0)
        || !d.lowest_excitation_hz.is_finite()
        || d.lowest_excitation_hz <= 0.0
        || d.lowest_excitation_hz > hi
    {
        return Err(binding(
            root,
            "claim",
            "finite nonnegative dynamics and positive excitation frequency no higher than claim band required",
        ));
    }
    let mut seen = BTreeSet::new();
    let mut limitations = Vec::new();
    for part in &input.parts {
        let mut covered = false;
        for m in claim.modules.iter().filter(|m| m.part == part.entity) {
            covered = true;
            let module = m.model.module();
            if !seen.insert((m.part, module)) {
                return Err(binding(m.part, module, "duplicate part/module"));
            }
            let row = claim
                .registry
                .row(&m.claim.filling, &m.claim.image, &m.claim.qoi)
                .filter(|row| row.gate != Gate::Refused)
                .ok_or_else(|| ValidityError::RegistryRefusal {
                    part: m.part,
                    claim: m.claim.clone(),
                })?;
            let expected_kind = match m.model {
                DescriptionModel::LumpedReed => PartKind::ValvePattern,
                DescriptionModel::OscillatingCylinderString => PartKind::String,
                _ => PartKind::Duct,
            };
            if part.kind != expected_kind {
                return Err(ValidityError::UnsupportedPart {
                    part: part.entity,
                    kind: part.kind,
                });
            }
            let mut states = vec![(None, part.clone())];
            for state in &input.configurations {
                let mut configured = part.clone();
                for (path, value) in &state.assignments {
                    if path.entity == part.entity {
                        let GeometrySource::Parametric { params, .. } = &mut configured.geometry
                        else {
                            return Err(ValidityError::UnsupportedPart {
                                part: part.entity,
                                kind: part.kind,
                            });
                        };
                        let (_, parameter) = params
                            .iter_mut()
                            .find(|(key, _)| *key == path.parameter)
                            .ok_or_else(|| missing(part, module, &path.parameter.name))?;
                        *parameter = *value;
                    }
                }
                states.push((Some(state.name.clone()), configured));
            }
            for (configuration, configured) in &states {
                let part = configured;
                let material_states = materials(part, m)?;
                let mut axes = vec![
                    policy(part, m, "frequency_hz", lo, HZ, None)?,
                    policy(part, m, "frequency_hz", hi, HZ, None)?,
                ];
                let statement = model_axes(part, m, claim, &mut axes)?;
                // Intersect declared policy with owner limits, even for axes
                // already evaluated above. Extra constrained axes are mandatory
                // description parameters, never implicitly unconstrained.
                for name in m.domain.bounds().keys() {
                    let mut values: Vec<_> = axes
                        .iter()
                        .filter(|a| &a.axis == name)
                        .map(|a| (a.value, a.unit, a.refinement))
                        .collect();
                    if values.is_empty() {
                        let (_, table) = params(part)?;
                        let (key, _) = table
                            .iter()
                            .find(|(k, _)| &k.name == name)
                            .ok_or_else(|| missing(part, module, name))?;
                        values.push((parameter(part, module, name, key.dims)?, key.dims, None));
                    }
                    for (value, unit, refinement) in values {
                        let a = policy(part, m, name, value, unit, refinement)?;
                        if !axes.contains(&a) {
                            axes.push(a);
                        }
                    }
                }
                limitations.push(ModelFormLimitation {
                    part: m.part,
                    module,
                    configuration: configuration.clone(),
                    statement,
                    claim: row.clone(),
                    axes,
                    material_states,
                });
            }
        }
        if !covered {
            return Err(ValidityError::UnsupportedPart {
                part: part.entity,
                kind: part.kind,
            });
        }
    }
    if input.parts.is_empty()
        || claim
            .modules
            .iter()
            .any(|m| !input.parts.iter().any(|p| p.entity == m.part))
    {
        return Err(binding(
            root,
            "claim",
            "nonempty exact part/model coverage required",
        ));
    }
    Ok(Admission {
        description,
        claim,
        limitations,
    })
}

fn model_axes(
    part: &PartDecl,
    m: &ModuleDomain<'_>,
    claim: &ClaimDomain<'_>,
    axes: &mut Vec<DomainAxis>,
) -> Result<&'static str, ValidityError> {
    let module = m.model.module();
    let gas = &claim.gas;
    let d = claim.dynamics;
    let omega = TAU * claim.band_hz.1;
    let positive = |name, value, unit| {
        axis(
            part,
            module,
            name,
            value,
            (f64::from_bits(1), f64::MAX),
            unit,
            None,
        )
    };
    match m.model {
        DescriptionModel::CylindricalDuct { loss } => {
            geometry(part, module, "cylinder.v1")?;
            let radius = parameter(part, module, "radius", LENGTH)?;
            let length = parameter(part, module, "length", LENGTH)?;
            Segment::Cylinder { radius, length }
                .validate()
                .map_err(|e| owner(part, module, e))?;
            // Rounded cutoff factors are the explicit §1.2 audit convention;
            // not a new numerical Bessel-root implementation.
            for (name, root) in [
                ("first_asymmetric_cutoff_hz", 1.841),
                ("first_axisymmetric_cutoff_hz", 3.832),
            ] {
                axes.push(axis(
                    part,
                    module,
                    name,
                    claim.band_hz.1,
                    (0.0, root * gas.sound_speed / (TAU * radius)),
                    HZ,
                    Some("M2.1-W3c / M2.4"),
                )?);
            }
            let wave_frequency = if loss == LossModel::WideTube {
                claim.band_hz.0
            } else {
                d.lowest_excitation_hz
            };
            if wave_frequency == 0.0 {
                axes.push(axis(
                    part,
                    module,
                    "shear_number",
                    0.0,
                    (fs_duct::MIN_SHEAR_NUMBER, f64::MAX),
                    Dims::NONE,
                    Some("M2.1-W5"),
                )?);
            }
            let wave = fs_duct::segment_wave(gas, radius, TAU * wave_frequency, loss)
                .map_err(|e| owner(part, module, e))?;
            if loss == LossModel::WideTube {
                axes.push(axis(
                    part,
                    module,
                    "shear_number",
                    wave.shear_number,
                    (fs_duct::MIN_SHEAR_NUMBER, f64::MAX),
                    Dims::NONE,
                    Some("M2.1-W5"),
                )?);
            }
            let evanescent = parameter(part, module, "evanescent_relative_error", Dims::NONE)?;
            axes.push(axis(
                part,
                module,
                "evanescent_relative_error",
                evanescent,
                (0.0, f64::MAX),
                Dims::NONE,
                Some("M2.4"),
            )?);
            axes.push(policy(
                part,
                m,
                "evanescent_relative_error",
                evanescent,
                Dims::NONE,
                Some("M2.4"),
            )?);
            // Inverse shock distance: beta*omega*p/(rho*c^3), beta=(gamma+1)/2.
            let ratio = length * ((gas.gamma + 1.0) * 0.5) * omega * d.max_mouth_pressure_pa
                / (gas.density * gas.sound_speed * gas.sound_speed * gas.sound_speed);
            axes.push(policy(
                part,
                m,
                "length_over_shock_distance",
                ratio,
                Dims::NONE,
                Some("M2.3"),
            )?);
            Ok(
                "1-D cylindrical linear duct; below transverse cutoff does not remove evanescent-field error. Declared error and weak-nonlinearity bounds are retained, not specimen qualification.",
            )
        }
        DescriptionModel::CompactToneHole => {
            geometry(part, module, "tone-hole.v1")?;
            let hole_radius = parameter(part, module, "hole_radius", LENGTH)?;
            let bore_radius = parameter(part, module, "bore_radius", LENGTH)?;
            let chimney_height = parameter(part, module, "chimney_height", LENGTH)?;
            let segment = Segment::ToneHole {
                hole_radius,
                bore_radius,
                chimney_height,
                state: fs_duct::HoleState::Open,
            };
            if let Err(e) = segment.validate() {
                if hole_radius >= bore_radius {
                    return Err(ValidityError::OutOfDomain {
                        part: part.entity,
                        axis: DomainAxis {
                            module,
                            axis: "hole_radius".into(),
                            value: hole_radius,
                            domain: (0.0, bore_radius),
                            upper_inclusive: false,
                            unit: LENGTH,
                            refinement: Some("M2.4"),
                        },
                    });
                }
                return Err(owner(part, module, e));
            }
            let edge_radius = parameter(part, module, "edge_radius", LENGTH)?;
            axes.push(positive("edge_radius", edge_radius, LENGTH)?);
            let displacement = d.max_mouth_pressure_pa
                / (gas.density * gas.sound_speed * TAU * d.lowest_excitation_hz);
            axes.push(policy(
                part,
                m,
                "edge_displacement_ratio",
                displacement / edge_radius,
                Dims::NONE,
                Some("M2.4 / M3.2"),
            )?);
            Ok(
                "Compact tone-hole T-junction; no resolved pad/finger/undercut or inter-hole exterior field. Edge-linearity envelope is declared, not a resolved flow certificate.",
            )
        }
        DescriptionModel::CompactRadiation { flange } => {
            geometry(part, module, "circular-opening.v1")?;
            let radius = parameter(part, module, "radius", LENGTH)?;
            axes.push(positive("radius", radius, LENGTH)?);
            let ka = omega / gas.sound_speed * radius;
            axes.push(axis(
                part,
                module,
                "ka",
                ka,
                (0.0, fs_duct::MAX_RADIATION_KA),
                Dims::NONE,
                Some("M2.5"),
            )?);
            // Keep the owner's own check in addition to a stricter declared box.
            fs_phs::compact_radiation_impedance(
                gas.density,
                gas.sound_speed,
                radius,
                omega,
                flange,
            )
            .map_err(|e| owner(part, module, e))?;
            axes.push(policy(part, m, "ka", ka, Dims::NONE, Some("M2.5"))?);
            Ok(
                "Compact low-ka resistive/reactive mouth fit; no full exterior directivity, bell field or arbitrary aperture shape.",
            )
        }
        DescriptionModel::LumpedReed => {
            geometry(part, module, "lumped-reed.v1")?;
            let p = |name, dims| parameter(part, module, name, dims);
            let reed = BeatingReed {
                rest_opening_m: p("rest_opening", LENGTH)?,
                width_m: p("width", LENGTH)?,
                closing_pressure_pa: p("closing_pressure", PRESSURE)?,
                blowing_pressure_pa: d.max_mouth_pressure_pa,
                attack_s: p("attack", Dims([0, 0, 1, 0, 0, 0]))?,
                mass_kg: p("mass", Dims([0, 1, 0, 0, 0, 0]))?,
                stiffness_n_m: p("stiffness", Dims([0, 1, -2, 0, 0, 0]))?,
                damping_ratio: p("damping_ratio", Dims::NONE)?,
            };
            if !crate::reed_bore::reed_parameters_valid(reed) {
                return Err(binding(
                    part.entity,
                    module,
                    "reed_parameters_valid refused",
                ));
            }
            axes.push(policy(
                part,
                m,
                "mouth_pressure_pa",
                d.max_mouth_pressure_pa,
                PRESSURE,
                Some("M3.1-W5 / M3.2"),
            )?);
            Ok(
                "Lumped reed/valve reduction; no resolved continuum, embouchure contact, anatomy or channel flow. Material resolution alone does not validate the reduction.",
            )
        }
        DescriptionModel::OscillatingCylinderString => {
            geometry(part, module, "circular-string.v1")?;
            let radius = parameter(part, module, "radius", LENGTH)?;
            axes.push(positive("radius", radius, LENGTH)?);
            // KC = Umax*T/D = 2*pi*displacement/(2*radius).
            axes.push(axis(
                part,
                module,
                "keulegan_carpenter",
                TAU * d.max_string_displacement_m / (2.0 * radius),
                (0.0, 1.0),
                Dims::NONE,
                Some("M3.2"),
            )?);
            axes.push(policy(
                part,
                m,
                "bow_force_n",
                d.max_bow_force_n,
                Dims([1, 1, -2, 0, 0, 0]),
                None,
            )?);
            axes.push(policy(
                part,
                m,
                "pluck_force_n",
                d.max_pluck_force_n,
                Dims([1, 1, -2, 0, 0, 0]),
                None,
            )?);
            Ok(
                "Oscillating-cylinder linear Stokes drag only at KC <= 1; no finite-amplitude separated flow or real-specimen qualification.",
            )
        }
    }
}

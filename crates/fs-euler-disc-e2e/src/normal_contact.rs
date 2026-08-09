//! Euler-disc adapter for the generic finite-patch normal-response laws.
//!
//! This is deliberately a narrow coordinate and admission adapter.  It maps a
//! bounded [`PatchKinematics`] record and caller-declared material/interface
//! inputs into `fs-contact` without selecting a material, fitting curvature,
//! or treating an event/barrier result as compliant contact.

use core::fmt;

use fs_contact::interface_binding::{BoundNormalContactLaw, BoundNormalContactModel};
use fs_contact::normal_patch::{
    ApplicabilityInput, ApplicabilityLimits, InputUncertainty, IntegrationLane,
    NormalPatchEmbedError, NormalPatchEmbedIdentity, NormalPatchEmbedRequest,
    NormalPatchEmbedState, NormalPatchEmbedTransition, NormalPatchGeometry, NormalPatchLaw,
    NormalPatchRequest,
};
use fs_mbd::Vec3;
use fs_tribo::{InputAuthority, InterfaceSystemRef};

use crate::patch_kinematics::{CurvatureMetadata, PatchContactStatus, PatchKinematics};

/// Stable identity of this coordinate-only normal-contact adapter.
pub const NORMAL_CONTACT_ADAPTER_ID: &str = "euler-disc/normal-contact-adapter-v1";

const CURVATURE_TOLERANCE: f64 = 256.0 * f64::EPSILON;

/// One caller-declared material and ordered interface input set.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalMaterialInterface {
    /// Stable material/interface-state binding identity.
    pub material_card_id: String,
    /// Stable normal-law identity retained in the generic receipt.
    pub model_id: String,
    /// Stable source-card identity retained in the generic receipt.
    pub source_id: String,
    /// Ordered interface/history/provenance data, forwarded without promotion.
    pub interface: InterfaceSystemRef,
    /// Caller-declared reduced modulus in Pa.
    pub reduced_modulus_pa: f64,
    /// Card-selected normal force-rate family. Geometry remains independently
    /// resolved from the actual patch and is never inferred from this field.
    pub rate_response: NormalRateResponse,
    /// Generic half-space, yield, rate, temperature, layer, and adhesion data.
    pub applicability: ApplicabilityInput,
    /// Explicit limits for the generic applicability ratios.
    pub limits: ApplicabilityLimits,
    /// Material/load input uncertainty; curvature uncertainty is merged in.
    pub uncertainty: InputUncertainty,
}

/// Geometry-owned settings needed to adapt an admitted material/model binding
/// to one normal-contact law. Bulk properties, ordered interface properties,
/// model choice, damping, rate scale, and applicability limits are
/// deliberately absent: the binding derives those from
/// [`BoundNormalContactModel`].
#[derive(Debug, Clone, PartialEq)]
pub struct NormalMaterialLawConfig {
    /// Local half-space depth used by the applicability check [m].
    pub half_space_depth_m: f64,
    /// Smallest relevant material-layer thickness [m].
    pub layer_thickness_m: f64,
    /// Propagated material/load uncertainty supplied by its evidence owner.
    pub uncertainty: InputUncertainty,
}

/// Build the Euler normal-contact material input from a shared, evidence-bound
/// ordered material state.
///
/// The temperature is the exact common `T` coordinate used to resolve both
/// material cards. The Hertz modulus and limiting yield stress are derived by
/// `fs-material`/`fs-contact`; this adapter cannot accept caller-retyped
/// duplicates of those values.
pub fn bind_normal_material_interface(
    model: &BoundNormalContactModel,
    config: NormalMaterialLawConfig,
) -> Result<NormalMaterialInterface, NormalMaterialBindingError> {
    let normal = model.normal_state();
    let elastic = normal.elastic();
    let temperature_k = elastic
        .state_coordinate("T")
        .ok_or(NormalMaterialBindingError::MissingTemperature)?;
    if !temperature_k.is_finite() || temperature_k <= 0.0 {
        return Err(NormalMaterialBindingError::InvalidTemperature { temperature_k });
    }
    for (value, field, strictly_positive) in [
        (config.half_space_depth_m, "half_space_depth_m", true),
        (config.layer_thickness_m, "layer_thickness_m", true),
    ] {
        if !value.is_finite()
            || (strictly_positive && value <= 0.0)
            || (!strictly_positive && value < 0.0)
        {
            return Err(NormalMaterialBindingError::InvalidConfig { field });
        }
    }
    let binding_id = model.identity();
    let model_card = model.model_card();
    let rate_response = match model.law() {
        BoundNormalContactLaw::ElasticHertz => NormalRateResponse::ElasticHertz,
        BoundNormalContactLaw::HuntCrossleyPoint {
            dissipation_s_per_m,
        } => NormalRateResponse::HuntCrossleyPoint {
            dissipation_s_per_m,
        },
    };
    Ok(NormalMaterialInterface {
        material_card_id: format!("fs-contact/normal-model-binding/{binding_id}"),
        model_id: format!(
            "{}:v{}:{}",
            model_card.law.0,
            model_card.law_version,
            model_card.content_hash()
        ),
        source_id: format!(
            "fs-matdb/constitutive-model-card/{}",
            model_card.content_hash()
        ),
        interface: elastic.interface().interface().clone(),
        reduced_modulus_pa: elastic.reduced_modulus_pa(),
        rate_response,
        applicability: ApplicabilityInput {
            half_space_depth_m: config.half_space_depth_m,
            layer_thickness_m: config.layer_thickness_m,
            yield_strength_pa: elastic.limiting_yield_stress_pa(),
            characteristic_rate_m_per_s: model.characteristic_rate_m_per_s(),
            temperature_k,
            adhesion_energy_j_per_m2: normal.adhesion_energy_j_per_m2(),
        },
        limits: model.limits(),
        uncertainty: config.uncertainty,
    })
}

/// Typed refusal from the ordered material-state to Euler normal-law bridge.
#[derive(Debug, Clone, PartialEq)]
pub enum NormalMaterialBindingError {
    /// The common material state point did not contain the canonical `T` axis.
    MissingTemperature,
    /// The retained absolute temperature was nonphysical.
    InvalidTemperature {
        /// Offered absolute temperature [K].
        temperature_k: f64,
    },
    /// A law-specific scalar or identity was invalid.
    InvalidConfig {
        /// Stable offending field name.
        field: &'static str,
    },
}

impl fmt::Display for NormalMaterialBindingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for NormalMaterialBindingError {}

/// The only local curvature shapes admitted by the current generic ladder.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EulerNormalGeometry {
    /// Both retained principal curvatures must agree within their retained
    /// geometry uncertainty; no caller-selected sphere fit is admitted.
    SpherePlane,
    /// One retained principal curvature must be flat within retained geometry
    /// uncertainty; the supplied line load is per unit axial length in N/m.
    CylinderPlane {
        /// Caller/outer-solver resolved line load in N/m.
        line_load_n_per_m: f64,
    },
    /// Both retained principal curvatures must be strictly positive *relative
    /// gap* curvatures, including the base contribution. The adapter cannot
    /// derive that pair from a disc-surface curvature alone. Its full ordered
    /// pair is passed to the elliptic Hertz law; no scalar effective-radius
    /// constitutive approximation is made.
    EllipticParaboloid,
}

/// Card-bound rate dependence of an admitted point normal response.
///
/// This is constitutive data, not a material-name dispatch. The actual local
/// curvature geometry remains an independent input to the normal law.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NormalRateResponse {
    /// Rate-independent elastic Hertz response.
    ElasticHertz,
    /// Passive Hunt--Crossley force factor for an admitted point-contact Hertz
    /// geometry. The coefficient has units s/m.
    HuntCrossleyPoint {
        /// Velocity-proportional indentation dissipation coefficient [s/m].
        dissipation_s_per_m: f64,
    },
}

/// Numerical regime owning the current normal-contact sample.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormalContactIntegrationRegime {
    /// Quasistatic/smooth fixed-branch solve; impact-candidate kinematics are a
    /// typed handoff rather than silently admitted.
    SmoothQuasistatic,
    /// Time-resolved compliant transient. High closure rates may be evaluated
    /// by the same finite-patch law, subject to its explicit rate, pressure,
    /// strain, temperature, and geometry applicability limits.
    CompliantTransient,
}

/// Identities owned by the Euler caller for one smooth fixed-branch sample.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NormalContactIdentity {
    /// Caller case/run identity.
    pub case_id: String,
    /// Must equal [`NORMAL_CONTACT_ADAPTER_ID`].
    pub adapter_id: String,
    /// Generic solver identity.
    pub solver_id: String,
    /// Generic contact identity.
    pub contact_id: String,
    /// Generic sample identity.
    pub sample_id: String,
}

/// Input for one Euler normal-contact evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct EulerNormalContactInput {
    /// Euler caller identities.
    pub identity: NormalContactIdentity,
    /// Bounded pre-constitutive geometry and relative motion.
    pub kinematics: PatchKinematics,
    /// Admitted-by-caller material/interface records.
    pub material: NormalMaterialInterface,
    /// Explicitly selected supported local curvature shape.
    pub geometry: EulerNormalGeometry,
    /// Numerical regime that owns this sample.
    pub integration_regime: NormalContactIntegrationRegime,
    /// Immutable generic exactly-once state from the prior accepted sample.
    pub state: NormalPatchEmbedState,
    /// Sample clock in seconds.
    pub time_s: f64,
    /// Solver iteration for the current smooth branch.
    pub iteration: u64,
    /// Accepted candidate interval duration in seconds.
    pub step_s: f64,
    /// A generic port is published only after the outer smooth iteration converges.
    pub converged: bool,
}

/// Retained principal-curvature decision; no effective-radius fit is hidden.
#[derive(Debug, Clone, PartialEq)]
pub struct CurvatureResolution {
    /// Geometry-owner curvature identity.
    pub curvature_identity: String,
    /// Geometry owner's authority, retained independently of material-card authority.
    pub authority: InputAuthority,
    /// First caller-supplied relative-gap principal curvature in 1/m.
    pub first_principal_m_inverse: f64,
    /// Second caller-supplied relative-gap principal curvature in 1/m.
    pub second_principal_m_inverse: f64,
    /// Retained absolute curvature uncertainty in 1/m.
    pub uncertainty_m_inverse: f64,
    /// Scalar radius retained for reporting and uncertainty normalization.
    /// The elliptic constitutive law does not consume this value.
    pub reporting_radius_m: f64,
}

/// Unit-typed elastic storage copied from a generic physical receipt.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NormalElasticStorage {
    /// Area/point contact reversible storage in J.
    PointJoules(f64),
    /// Line contact reversible storage normalized by axial length in J/m.
    LineJoulesPerMetre(f64),
}

/// Unit-typed irreversible work copied from a generic physical receipt.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NormalDissipation {
    /// Point-contact irreversible work and instantaneous dissipated power.
    Point { work_j: f64, power_w: f64 },
    /// Line-contact irreversible work and power, each normalized by length.
    Line {
        work_j_per_m: f64,
        power_w_per_m: f64,
    },
}

/// An accepted, active finite-patch normal-contact mapping.
#[derive(Debug, Clone, PartialEq)]
pub struct ActiveNormalContact {
    /// Exact generic receipt, generic port, applicability, uncertainty, and successor state.
    pub generic: NormalPatchEmbedTransition,
    /// Retained material-card identity; this adapter neither admits nor rewrites it.
    pub material_card_id: String,
    /// Retained caller source-card identity.
    pub material_source_id: String,
    /// Disc application point in world metres.
    pub application_point_world_m: Vec3,
    /// Disc-centre-of-mass-relative application arm in world metres.
    pub application_arm_world_m: Vec3,
    /// Both source principal curvatures and the explicit supported reduction.
    pub curvature: CurvatureResolution,
    /// Reversible elastic storage with point/line units preserved.
    pub elastic_storage: NormalElasticStorage,
    /// Irreversible work/power with point/line units preserved.
    pub dissipation: NormalDissipation,
}

/// Explicit non-active or active outcome.  Separation never creates a zero-force port.
#[derive(Debug, Clone, PartialEq)]
pub enum EulerNormalContactOutcome {
    /// Contact is certainly separated; the state is unmodified and no law was evaluated.
    InactiveSeparated {
        /// Positive-opening gap retained from kinematics in metres.
        gap_m: f64,
        /// Positive means opening, retained in m/s.
        normal_relative_velocity_m_per_s: f64,
        /// Unchanged generic state; no work key was consumed.
        state: NormalPatchEmbedState,
    },
    /// A smooth fixed-branch physical compliant response was admitted.
    Active(ActiveNormalContact),
}

/// Refusal surface for Euler-specific mapping before a port is published.
#[derive(Debug, Clone, PartialEq)]
pub enum NormalContactError {
    /// A required identity was blank or did not bind this adapter version.
    InvalidIdentity { field: &'static str },
    /// A scalar/vector was non-finite or outside this adapter's domain.
    InvalidInput { field: &'static str },
    /// Unknown and impact-candidate kinematics must be resolved by their owner.
    UnavailableKinematics { status: PatchContactStatus },
    /// Curvature was withheld by its geometry owner.
    CurvatureUnavailable,
    /// The two supplied principal curvatures do not meet the declared sphere condition.
    SphereCurvatureMismatch {
        first_m_inverse: f64,
        second_m_inverse: f64,
        tolerance_m_inverse: f64,
    },
    /// The two supplied principal curvatures do not meet the declared cylinder condition.
    CylinderCurvatureMismatch {
        first_m_inverse: f64,
        second_m_inverse: f64,
        flatness_tolerance_m_inverse: f64,
    },
    /// A two-curvature/high-ellipticity patch needs an independently authorized law.
    UnsupportedTwoCurvature {
        first_m_inverse: f64,
        second_m_inverse: f64,
    },
    /// Hunt--Crossley is not a line-contact law in the generic ladder.
    DissipativeLineUnsupported,
    /// The generic finite-patch law or transactional embedding refused the request.
    GenericRefusal(NormalPatchEmbedError),
}

impl fmt::Display for NormalContactError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for NormalContactError {}

impl From<NormalPatchEmbedError> for NormalContactError {
    fn from(value: NormalPatchEmbedError) -> Self {
        Self::GenericRefusal(value)
    }
}

/// Maps one bounded Euler patch into the generic finite-patch normal law.
///
/// The local normal is directed base-to-disc, matching the action on the disc.
/// `PatchKinematics::normal_relative_velocity_m_per_s` is positive while the
/// gap opens, so the generic indentation rate is its negation.
pub fn evaluate_normal_contact(
    input: &EulerNormalContactInput,
) -> Result<EulerNormalContactOutcome, NormalContactError> {
    validate_identity(&input.identity)?;
    validate_kinematics(&input.kinematics)?;
    let gap_m = normal_gap_m(&input.kinematics)?;

    match input.kinematics.status {
        PatchContactStatus::Separated => {
            return Ok(EulerNormalContactOutcome::InactiveSeparated {
                gap_m,
                normal_relative_velocity_m_per_s: input.kinematics.normal_relative_velocity_m_per_s,
                state: input.state.clone(),
            });
        }
        PatchContactStatus::Unknown => {
            return Err(NormalContactError::UnavailableKinematics {
                status: input.kinematics.status,
            });
        }
        PatchContactStatus::ImpactCandidate | PatchContactStatus::Receding
            if input.integration_regime != NormalContactIntegrationRegime::CompliantTransient =>
        {
            return Err(NormalContactError::UnavailableKinematics {
                status: input.kinematics.status,
            });
        }
        PatchContactStatus::ImpactCandidate | PatchContactStatus::Receding => {}
        PatchContactStatus::Approaching
        | PatchContactStatus::Touching
        | PatchContactStatus::Grazing => {}
    }

    validate_sample(input)?;
    validate_material(&input.material)?;
    let curvature = resolve_curvature(&input.kinematics, input.geometry)?;
    let normal = input.kinematics.tangent_basis.normal_world;
    let approach_rate_m_per_s = -input.kinematics.normal_relative_velocity_m_per_s;
    let (law, geometry, line_load_n_per_m) =
        law_from_input(&input.material, input.geometry, &curvature)?;
    let uncertainty = merged_uncertainty(input.material.uncertainty, &curvature)?;
    let request = NormalPatchRequest {
        identity: fs_contact::normal_patch::NormalPatchIdentity {
            model_id: input.material.model_id.clone(),
            source_id: input.material.source_id.clone(),
            state_id: input.identity.case_id.clone(),
        },
        interface: input.material.interface.clone(),
        law,
        geometry,
        indentation_m: (-gap_m).max(0.0),
        indentation_rate_m_per_s: approach_rate_m_per_s,
        step_s: input.step_s,
        line_load_n_per_m,
        applicability: input.material.applicability,
        limits: input.material.limits,
        uncertainty,
    };
    let generic = NormalPatchEmbedRequest {
        identity: NormalPatchEmbedIdentity {
            solver_id: input.identity.solver_id.clone(),
            contact_id: input.identity.contact_id.clone(),
            feature_id: input.kinematics.patch.patch_identity.as_str().to_owned(),
            sample_id: input.identity.sample_id.clone(),
        },
        lane: IntegrationLane::SmoothFixed,
        converged: input.converged,
        kinematics: fs_contact::normal_patch::NormalPatchKinematics {
            declared_gap_m: gap_m,
            approach_m: (-gap_m).max(0.0),
            approach_rate_m_per_s,
            time_s: input.time_s,
            step_s: input.step_s,
            iteration: input.iteration,
            normal: [normal.x, normal.y, normal.z],
            moment_arm_m: [
                input.kinematics.disc_point.arm_world.x,
                input.kinematics.disc_point.arm_world.y,
                input.kinematics.disc_point.arm_world.z,
            ],
        },
        law_request: request,
    }
    .evaluate(&input.state)?;
    let (elastic_storage, dissipation) = storage_and_dissipation(&generic);
    Ok(EulerNormalContactOutcome::Active(ActiveNormalContact {
        application_point_world_m: input.kinematics.disc_point.point_world,
        application_arm_world_m: input.kinematics.disc_point.arm_world,
        material_card_id: input.material.material_card_id.clone(),
        material_source_id: input.material.source_id.clone(),
        curvature,
        elastic_storage,
        dissipation,
        generic,
    }))
}

fn validate_identity(identity: &NormalContactIdentity) -> Result<(), NormalContactError> {
    for (value, field) in [
        (identity.case_id.as_str(), "case_id"),
        (identity.solver_id.as_str(), "solver_id"),
        (identity.contact_id.as_str(), "contact_id"),
        (identity.sample_id.as_str(), "sample_id"),
    ] {
        if value.trim().is_empty() {
            return Err(NormalContactError::InvalidIdentity { field });
        }
    }
    if identity.adapter_id != NORMAL_CONTACT_ADAPTER_ID {
        return Err(NormalContactError::InvalidIdentity {
            field: "adapter_id",
        });
    }
    Ok(())
}

fn validate_material(material: &NormalMaterialInterface) -> Result<(), NormalContactError> {
    for (value, field) in [
        (material.material_card_id.as_str(), "material_card_id"),
        (material.model_id.as_str(), "model_id"),
        (material.source_id.as_str(), "source_id"),
    ] {
        if value.trim().is_empty() {
            return Err(NormalContactError::InvalidIdentity { field });
        }
    }
    if !material.reduced_modulus_pa.is_finite() || material.reduced_modulus_pa <= 0.0 {
        return Err(NormalContactError::InvalidInput {
            field: "reduced_modulus_pa",
        });
    }
    if let NormalRateResponse::HuntCrossleyPoint {
        dissipation_s_per_m,
    } = material.rate_response
        && (!dissipation_s_per_m.is_finite() || dissipation_s_per_m < 0.0)
    {
        return Err(NormalContactError::InvalidInput {
            field: "hunt_crossley_dissipation_s_per_m",
        });
    }
    Ok(())
}

fn validate_sample(input: &EulerNormalContactInput) -> Result<(), NormalContactError> {
    if !input.time_s.is_finite() || input.time_s < 0.0 {
        return Err(NormalContactError::InvalidInput { field: "time_s" });
    }
    if !input.step_s.is_finite() || input.step_s <= 0.0 {
        return Err(NormalContactError::InvalidInput { field: "step_s" });
    }
    validate_kinematics(&input.kinematics)
}

fn validate_kinematics(kinematics: &PatchKinematics) -> Result<(), NormalContactError> {
    if !kinematics.normal_relative_velocity_m_per_s.is_finite()
        || !kinematics.tangent_basis.normal_world.is_finite()
        || !kinematics.disc_point.point_world.is_finite()
        || !kinematics.disc_point.arm_world.is_finite()
    {
        return Err(NormalContactError::InvalidInput {
            field: "patch kinematics",
        });
    }
    Ok(())
}

fn normal_gap_m(kinematics: &PatchKinematics) -> Result<f64, NormalContactError> {
    let gap_m = kinematics
        .disc_point
        .point_world
        .sub(kinematics.base_point.point_world)
        .dot(kinematics.tangent_basis.normal_world);
    if gap_m.is_finite() {
        Ok(gap_m)
    } else {
        Err(NormalContactError::InvalidInput {
            field: "patch gap reconstruction",
        })
    }
}

fn resolve_curvature(
    kinematics: &PatchKinematics,
    geometry: EulerNormalGeometry,
) -> Result<CurvatureResolution, NormalContactError> {
    // `PatchKinematics` is a geometry-owner boundary. For the elliptic law,
    // its retained values must already be principal curvatures of the local
    // *relative gap* (disc plus base); this adapter has neither a second
    // surface chart nor authority to manufacture the missing base curvature.
    let CurvatureMetadata::Known {
        curvature_identity,
        authority,
        first_principal_m_inverse,
        second_principal_m_inverse,
        uncertainty_m_inverse,
    } = &kinematics.patch.curvature
    else {
        return Err(NormalContactError::CurvatureUnavailable);
    };
    let first = *first_principal_m_inverse;
    let second = *second_principal_m_inverse;
    let uncertainty = *uncertainty_m_inverse;
    if !first.is_finite() || !second.is_finite() || !uncertainty.is_finite() || uncertainty < 0.0 {
        return Err(NormalContactError::InvalidInput {
            field: "principal curvature",
        });
    }
    let reporting_radius_m = match geometry {
        EulerNormalGeometry::SpherePlane => {
            if first <= 0.0 || second <= 0.0 {
                return Err(NormalContactError::UnsupportedTwoCurvature {
                    first_m_inverse: first,
                    second_m_inverse: second,
                });
            }
            let tolerance = uncertainty.max(CURVATURE_TOLERANCE);
            if (first - second).abs() > tolerance {
                return Err(NormalContactError::SphereCurvatureMismatch {
                    first_m_inverse: first,
                    second_m_inverse: second,
                    tolerance_m_inverse: tolerance,
                });
            }
            2.0 / (first + second)
        }
        EulerNormalGeometry::CylinderPlane { .. } => {
            let flat = uncertainty.max(CURVATURE_TOLERANCE);
            let curved = if first.abs() <= flat && second > 0.0 {
                second
            } else if second.abs() <= flat && first > 0.0 {
                first
            } else {
                return Err(NormalContactError::CylinderCurvatureMismatch {
                    first_m_inverse: first,
                    second_m_inverse: second,
                    flatness_tolerance_m_inverse: flat,
                });
            };
            1.0 / curved
        }
        EulerNormalGeometry::EllipticParaboloid => {
            if first <= 0.0 || second <= 0.0 {
                return Err(NormalContactError::UnsupportedTwoCurvature {
                    first_m_inverse: first,
                    second_m_inverse: second,
                });
            }
            1.0 / (first * second).sqrt()
        }
    };
    if !reporting_radius_m.is_finite() || reporting_radius_m <= 0.0 {
        return Err(NormalContactError::InvalidInput {
            field: "reporting_radius_m",
        });
    }
    Ok(CurvatureResolution {
        curvature_identity: curvature_identity.as_str().to_owned(),
        authority: *authority,
        first_principal_m_inverse: first,
        second_principal_m_inverse: second,
        uncertainty_m_inverse: uncertainty,
        reporting_radius_m,
    })
}

fn law_from_input(
    material: &NormalMaterialInterface,
    geometry: EulerNormalGeometry,
    curvature: &CurvatureResolution,
) -> Result<(NormalPatchLaw, NormalPatchGeometry, f64), NormalContactError> {
    match geometry {
        EulerNormalGeometry::SpherePlane => Ok((
            match material.rate_response {
                NormalRateResponse::HuntCrossleyPoint {
                    dissipation_s_per_m,
                } => NormalPatchLaw::HuntCrossleySphere {
                    effective_radius_m: curvature.reporting_radius_m,
                    reduced_modulus_pa: material.reduced_modulus_pa,
                    dissipation_s_per_m,
                },
                NormalRateResponse::ElasticHertz => NormalPatchLaw::HertzSpherePlane {
                    effective_radius_m: curvature.reporting_radius_m,
                    reduced_modulus_pa: material.reduced_modulus_pa,
                },
            },
            NormalPatchGeometry::SpherePlane,
            0.0,
        )),
        EulerNormalGeometry::CylinderPlane {
            line_load_n_per_m, ..
        } => {
            if matches!(
                material.rate_response,
                NormalRateResponse::HuntCrossleyPoint { .. }
            ) {
                return Err(NormalContactError::DissipativeLineUnsupported);
            }
            if !line_load_n_per_m.is_finite() || line_load_n_per_m < 0.0 {
                return Err(NormalContactError::InvalidInput {
                    field: "line_load_n_per_m",
                });
            }
            Ok((
                NormalPatchLaw::HertzCylinderPlane {
                    effective_radius_m: curvature.reporting_radius_m,
                    reduced_modulus_pa: material.reduced_modulus_pa,
                },
                NormalPatchGeometry::CylinderPlane,
                line_load_n_per_m,
            ))
        }
        EulerNormalGeometry::EllipticParaboloid => {
            let maximum_principal_curvature_m_inverse = curvature
                .first_principal_m_inverse
                .max(curvature.second_principal_m_inverse);
            let minimum_principal_curvature_m_inverse = curvature
                .first_principal_m_inverse
                .min(curvature.second_principal_m_inverse);
            let law = match material.rate_response {
                NormalRateResponse::ElasticHertz => NormalPatchLaw::HertzEllipticParaboloid {
                    maximum_principal_curvature_m_inverse,
                    minimum_principal_curvature_m_inverse,
                    reduced_modulus_pa: material.reduced_modulus_pa,
                },
                NormalRateResponse::HuntCrossleyPoint {
                    dissipation_s_per_m,
                } => NormalPatchLaw::HuntCrossleyEllipticParaboloid {
                    maximum_principal_curvature_m_inverse,
                    minimum_principal_curvature_m_inverse,
                    reduced_modulus_pa: material.reduced_modulus_pa,
                    dissipation_s_per_m,
                },
            };
            Ok((law, NormalPatchGeometry::EllipticParaboloid, 0.0))
        }
    }
}

fn merged_uncertainty(
    material: InputUncertainty,
    curvature: &CurvatureResolution,
) -> Result<InputUncertainty, NormalContactError> {
    let positive_minimum = curvature
        .first_principal_m_inverse
        .min(curvature.second_principal_m_inverse);
    let curvature_scale = if positive_minimum > 0.0 {
        positive_minimum
    } else {
        (1.0 / curvature.reporting_radius_m).max(f64::MIN_POSITIVE)
    };
    let curvature_relative = curvature.uncertainty_m_inverse / curvature_scale;
    if !curvature_relative.is_finite() {
        return Err(NormalContactError::InvalidInput {
            field: "curvature uncertainty",
        });
    }
    Ok(InputUncertainty {
        radius_relative: material.radius_relative.max(curvature_relative),
        modulus_relative: material.modulus_relative,
        load_relative: material.load_relative,
    })
}

fn storage_and_dissipation(
    transition: &NormalPatchEmbedTransition,
) -> (NormalElasticStorage, NormalDissipation) {
    match &transition.receipt {
        fs_contact::normal_patch::NormalPatchReceipt::Point(receipt) => (
            NormalElasticStorage::PointJoules(receipt.reversible_energy_j),
            NormalDissipation::Point {
                work_j: receipt.irreversible_work_j,
                power_w: receipt.dissipated_power_w,
            },
        ),
        fs_contact::normal_patch::NormalPatchReceipt::Line(receipt) => (
            NormalElasticStorage::LineJoulesPerMetre(receipt.reversible_energy_j_per_m),
            NormalDissipation::Line {
                work_j_per_m: receipt.irreversible_work_j_per_m,
                power_w_per_m: receipt.dissipated_power_w_per_m,
            },
        ),
    }
}

#[cfg(test)]
mod material_binding_tests {
    use std::collections::BTreeMap;

    use fs_contact::interface_binding::{
        ADHESION_ENERGY_DIMS, ADHESION_ENERGY_PROPERTY, NORMAL_HERTZ_LAW_ID,
        NormalContactModelSelection, bind_dry_interface_system_card,
        bind_isotropic_elastic_interface, bind_normal_contact_model, bind_normal_interface_state,
    };
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterfaceSystemCard,
        InterpolationPolicy, LawId, LawParameter, MaterialCard, MaterialStateId, PropertyClaim,
        PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
        UncertaintyModel,
    };
    use fs_material::state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_interface_state_point, resolve_isotropic_solid_state_point,
    };
    use fs_qty::{Density, Dims, Pressure};

    use super::*;

    fn state_id(chemistry: &str) -> MaterialStateId {
        MaterialStateId {
            chemistry: chemistry.to_owned(),
            phase: "solid".to_owned(),
            process: "synthetic-state-series".to_owned(),
            revision: 0,
        }
    }

    fn material_card(
        chemistry: &str,
        density: f64,
        young: f64,
        poisson: f64,
        yield_stress: f64,
    ) -> MaterialCard {
        let mut claims = ClaimSet::new();
        for (name, dims, value) in [
            ("density", Density::DIMS, density),
            ("young_modulus", Pressure::DIMS, young),
            ("poisson_ratio", Dims::NONE, poisson),
            ("yield_stress", Pressure::DIMS, yield_stress),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Scalar { value, dims },
                    validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                    uncertainty: UncertaintyModel::RelativeHalfWidth {
                        fraction: 0.01,
                        confidence: 0.95,
                    },
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: Vec::new(),
                    provenance: Provenance {
                        source: format!("synthetic {chemistry} {name}"),
                        license: "CC0-1.0".to_owned(),
                        artifact: None,
                    },
                })
                .expect("property claim");
        }
        MaterialCard::assemble(state_id(chemistry), claims, Vec::new()).expect("material card")
    }

    fn normal_model_card() -> ConstitutiveModelCard {
        let mut parameters = BTreeMap::new();
        for (name, dims, value) in [
            ("characteristic-rate", Dims([1, 0, -1, 0, 0, 0]), 1.0),
            ("max-patch-to-radius", Dims::NONE, 0.1),
            ("max-strain", Dims::NONE, 0.01),
            ("max-patch-to-depth", Dims::NONE, 0.1),
            ("max-patch-to-layer", Dims::NONE, 0.1),
            ("max-pressure-to-yield", Dims::NONE, 0.2),
            ("max-rate-ratio", Dims::NONE, 0.1),
        ] {
            parameters.insert(name.to_owned(), LawParameter { value, dims });
        }
        ConstitutiveModelCard {
            law: LawId(NORMAL_HERTZ_LAW_ID.to_owned()),
            law_version: 1,
            parameters,
            state_schema_version: 1,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
            sources: Vec::new(),
            provenance: Provenance {
                source: "synthetic Hertz applicability card".to_owned(),
                license: "CC0-1.0".to_owned(),
                artifact: None,
            },
        }
    }

    #[test]
    fn g0_normal_material_binding_uses_resolved_cards_and_common_temperature() {
        let copper_card = material_card("copper-c110", 8960.0, 117.0e9, 0.34, 70.0e6);
        let glass_card = material_card("soda-lime-glass", 2500.0, 72.0e9, 0.22, 1.0e9);
        let point = QueryPoint::new().with("T", 293.15).expect("state point");
        let copper = resolve_isotropic_solid_state_point(
            &copper_card,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .expect("copper state");
        let glass = resolve_isotropic_solid_state_point(
            &glass_card,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .expect("glass state");
        let mut interface_claims = ClaimSet::new();
        interface_claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(ADHESION_ENERGY_PROPERTY, ADHESION_ENERGY_DIMS),
                value: PropertyValue::Scalar {
                    value: 0.0,
                    dims: ADHESION_ENERGY_DIMS,
                },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                uncertainty: UncertaintyModel::HalfWidth {
                    half_width: 0.0,
                    confidence: 0.95,
                },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: Vec::new(),
                provenance: Provenance {
                    source: "synthetic explicitly nonadhesive interface".to_owned(),
                    license: "CC0-1.0".to_owned(),
                    artifact: None,
                },
            })
            .expect("adhesion claim");
        let interface_card = InterfaceSystemCard::assemble(
            SurfaceSpec {
                material: state_id("copper-c110"),
                texture_frame: "disc-edge/profile-17".to_owned(),
            },
            SurfaceSpec {
                material: state_id("soda-lime-glass"),
                texture_frame: "base-track/profile-4".to_owned(),
            },
            SystemContext {
                medium: "dry".to_owned(),
                third_body: None,
                environment: "air-293K".to_owned(),
                history: "cleaned".to_owned(),
            },
            interface_claims,
            vec![normal_model_card()],
        )
        .expect("interface card");
        let dry = bind_dry_interface_system_card(&interface_card, InputAuthority::SyntheticFixture)
            .expect("dry interface");
        let elastic = bind_isotropic_elastic_interface(&dry, &copper, &glass)
            .expect("ordered elastic interface");
        let interface_state = resolve_interface_state_point(
            &interface_card,
            &point,
            &[ScalarPropertyRequirement::try_new(
                ADHESION_ENERGY_PROPERTY,
                ADHESION_ENERGY_DIMS,
                ScalarAdmissibility::NonNegative,
            )
            .expect("normal interface requirement")],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .expect("interface property state");
        let normal = bind_normal_interface_state(&elastic, &interface_state)
            .expect("normal interface state");
        let model =
            bind_normal_contact_model(&normal, NormalContactModelSelection::SingleSupported)
                .expect("normal model binding");
        let material = bind_normal_material_interface(
            &model,
            NormalMaterialLawConfig {
                half_space_depth_m: 0.02,
                layer_thickness_m: 0.01,
                uncertainty: InputUncertainty {
                    radius_relative: 0.001,
                    modulus_relative: 0.02,
                    load_relative: 0.03,
                },
            },
        )
        .expect("normal material binding");
        assert_eq!(
            material.reduced_modulus_pa.to_bits(),
            elastic.reduced_modulus_pa().to_bits()
        );
        assert_eq!(material.applicability.yield_strength_pa, 70.0e6);
        assert_eq!(material.applicability.temperature_k, 293.15);
        assert_eq!(material.applicability.characteristic_rate_m_per_s, 1.0);
        assert_eq!(material.limits.max_pressure_to_yield, 0.2);
        assert!(material.model_id.starts_with(NORMAL_HERTZ_LAW_ID));
        assert_eq!(
            material.interface.ordered_system_id(),
            dry.interface().ordered_system_id()
        );
    }
}

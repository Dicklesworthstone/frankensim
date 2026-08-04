//! Euler-disc adapter for the generic finite-patch normal-response laws.
//!
//! This is deliberately a narrow coordinate and admission adapter.  It maps a
//! bounded [`PatchKinematics`] record and caller-declared material/interface
//! inputs into `fs-contact` without selecting a material, fitting curvature,
//! or treating an event/barrier result as compliant contact.

use core::fmt;

use fs_contact::normal_patch::{
    ApplicabilityInput, ApplicabilityLimits, InputUncertainty, IntegrationLane,
    NormalPatchEmbedError, NormalPatchEmbedIdentity, NormalPatchEmbedRequest,
    NormalPatchEmbedState, NormalPatchEmbedTransition, NormalPatchGeometry, NormalPatchLaw,
    NormalPatchRequest,
};
use fs_mbd::Vec3;
use fs_tribo::InterfaceSystemRef;

use crate::patch_kinematics::{CurvatureMetadata, PatchContactStatus, PatchKinematics};

/// Stable identity of this coordinate-only normal-contact adapter.
pub const NORMAL_CONTACT_ADAPTER_ID: &str = "euler-disc/normal-contact-adapter-v1";

const CURVATURE_TOLERANCE: f64 = 256.0 * f64::EPSILON;

/// One caller-declared material and ordered interface input set.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalMaterialInterface {
    /// Stable material-card identity; this adapter does not admit the card.
    pub material_card_id: String,
    /// Stable normal-law identity retained in the generic receipt.
    pub model_id: String,
    /// Stable source-card identity retained in the generic receipt.
    pub source_id: String,
    /// Ordered interface/history/provenance data, forwarded without promotion.
    pub interface: InterfaceSystemRef,
    /// Caller-declared reduced modulus in Pa.
    pub reduced_modulus_pa: f64,
    /// Optional Hunt--Crossley dissipation coefficient in s/m.  It is valid
    /// only for the admitted sphere/plane rung.
    pub hunt_crossley_dissipation_s_per_m: Option<f64>,
    /// Generic half-space, yield, rate, temperature, layer, and adhesion data.
    pub applicability: ApplicabilityInput,
    /// Explicit limits for the generic applicability ratios.
    pub limits: ApplicabilityLimits,
    /// Material/load input uncertainty; curvature uncertainty is merged in.
    pub uncertainty: InputUncertainty,
}

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
    /// The current generic ladder has no dissipative elliptic Hertz variant.
    DissipativeEllipticUnsupported,
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
        PatchContactStatus::Unknown | PatchContactStatus::ImpactCandidate => {
            return Err(NormalContactError::UnavailableKinematics {
                status: input.kinematics.status,
            });
        }
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
    if let Some(value) = material.hunt_crossley_dissipation_s_per_m
        && (!value.is_finite() || value < 0.0)
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
            match material.hunt_crossley_dissipation_s_per_m {
                Some(dissipation_s_per_m) => NormalPatchLaw::HuntCrossleySphere {
                    effective_radius_m: curvature.reporting_radius_m,
                    reduced_modulus_pa: material.reduced_modulus_pa,
                    dissipation_s_per_m,
                },
                None => NormalPatchLaw::HertzSpherePlane {
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
            if material.hunt_crossley_dissipation_s_per_m.is_some() {
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
            if material.hunt_crossley_dissipation_s_per_m.is_some() {
                return Err(NormalContactError::DissipativeEllipticUnsupported);
            }
            let maximum_principal_curvature_m_inverse = curvature
                .first_principal_m_inverse
                .max(curvature.second_principal_m_inverse);
            let minimum_principal_curvature_m_inverse = curvature
                .first_principal_m_inverse
                .min(curvature.second_principal_m_inverse);
            Ok((
                NormalPatchLaw::HertzEllipticParaboloid {
                    maximum_principal_curvature_m_inverse,
                    minimum_principal_curvature_m_inverse,
                    reduced_modulus_pa: material.reduced_modulus_pa,
                },
                NormalPatchGeometry::EllipticParaboloid,
                0.0,
            ))
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

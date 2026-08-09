//! Analytic Hertz and Hunt--Crossley normal-patch laws in SI units.
//!
//! The physical receipt type is intentionally distinct from the nonphysical
//! constraint/barrier marker: an IPC/barrier result cannot be passed where a
//! compliance response is required.

use core::fmt;

use fs_blake3::{ContentHash, hash_domain};
use fs_tribo::{InputAuthority, InterfaceSystemRef};

/// SI dimensions for point-contact receipts.
pub const POINT_SI_UNITS: &str = "SI:m,s,N,N/m,Pa,J,W";
/// SI dimensions for line-contact receipts, normalized per unit cylinder length.
pub const LINE_SI_UNITS: &str = "SI:m,s,N/m,Pa,J/m,W/m";
const DOMAIN: &str = "org.frankensim.fs-contact.normal-patch.v1";

/// Explicit identity for a reproducible law request.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NormalPatchIdentity {
    /// Stable law/model identifier.
    pub model_id: String,
    /// Caller-provided source identifier.
    pub source_id: String,
    /// State/history identity.
    pub state_id: String,
}

/// Analytic normal law rung.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NormalPatchLaw {
    /// Elastic Hertz sphere against a half space.
    HertzSpherePlane {
        effective_radius_m: f64,
        reduced_modulus_pa: f64,
    },
    /// Elastic Hertz cylinder against a half space, per unit cylinder length.
    HertzCylinderPlane {
        effective_radius_m: f64,
        reduced_modulus_pa: f64,
    },
    /// Passive Hunt--Crossley augmentation of the Hertz sphere rung.
    HuntCrossleySphere {
        effective_radius_m: f64,
        reduced_modulus_pa: f64,
        dissipation_s_per_m: f64,
    },
    /// Elastic Hertz contact for the local separation paraboloid
    /// `0.5 * (k_max * x^2 + k_min * y^2)` against a half space.
    ///
    /// The curvatures are those of the *relative local gap*, not an inferred
    /// radius of either body. `maximum_principal_curvature_m_inverse` must be
    /// at least `minimum_principal_curvature_m_inverse`, and both must be
    /// strictly positive. The resulting patch is generally elliptical.
    /// Caller-declared uncertainty is retained in the receipt but this rung
    /// does not claim a propagated uncertainty bound.
    HertzEllipticParaboloid {
        maximum_principal_curvature_m_inverse: f64,
        minimum_principal_curvature_m_inverse: f64,
        reduced_modulus_pa: f64,
    },
    /// Passive Hunt--Crossley augmentation of the elliptic Hertz rung.
    ///
    /// This retains the true two-principal-curvature Hertz coefficient and
    /// multiplies its elastic force by `1 + c * indentation_rate`. It is a
    /// compliant transient point-contact law, not a coefficient-of-restitution
    /// shortcut. The request applicability envelope remains responsible for
    /// refusing rates, pressures, strains, or temperatures outside the model
    /// card's evidence.
    HuntCrossleyEllipticParaboloid {
        maximum_principal_curvature_m_inverse: f64,
        minimum_principal_curvature_m_inverse: f64,
        reduced_modulus_pa: f64,
        dissipation_s_per_m: f64,
    },
}

/// Declared local geometry. The analytic rungs deliberately do not approximate
/// toroidal or highly elliptical patches.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormalPatchGeometry {
    /// Axisymmetric sphere against a plane half space.
    SpherePlane,
    /// Cylinder against a plane half space, normalized per unit axial length.
    CylinderPlane,
    /// A locally quadratic, two-positive-principal-curvature separation.
    ///
    /// This is intentionally narrower than a generic toroidal body: callers
    /// must establish the local paraboloid and supply the relative curvatures.
    EllipticParaboloid,
    /// A two-curvature or strongly elliptical patch requiring another law.
    ToroidalOrHighlyElliptical,
}

/// Inputs governing the validity envelope. All ratios are dimensionless.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ApplicabilityInput {
    pub half_space_depth_m: f64,
    pub layer_thickness_m: f64,
    pub yield_strength_pa: f64,
    pub characteristic_rate_m_per_s: f64,
    pub temperature_k: f64,
    pub adhesion_energy_j_per_m2: f64,
}

/// Declared refusal limits for the small-strain, nonadhesive half-space law.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ApplicabilityLimits {
    pub max_patch_to_radius: f64,
    pub max_strain: f64,
    pub max_patch_to_depth: f64,
    pub max_patch_to_layer: f64,
    pub max_pressure_to_yield: f64,
    pub max_rate_ratio: f64,
    pub min_temperature_k: f64,
    pub max_temperature_k: f64,
}

/// Caller-declared relative input uncertainty; preserved without promotion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InputUncertainty {
    pub radius_relative: f64,
    pub modulus_relative: f64,
    pub load_relative: f64,
}

/// One finite, bounded request. `line_load_n_per_m` applies only to cylinder contact.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalPatchRequest {
    pub identity: NormalPatchIdentity,
    pub interface: InterfaceSystemRef,
    pub law: NormalPatchLaw,
    pub geometry: NormalPatchGeometry,
    pub indentation_m: f64,
    pub indentation_rate_m_per_s: f64,
    pub step_s: f64,
    pub line_load_n_per_m: f64,
    pub applicability: ApplicabilityInput,
    pub limits: ApplicabilityLimits,
    pub uncertainty: InputUncertainty,
}

/// Derived ratios retained in a successful physical receipt.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ApplicabilityRatios {
    pub patch_to_radius: f64,
    pub strain: f64,
    pub patch_to_depth: f64,
    pub patch_to_layer: f64,
    pub pressure_to_yield: f64,
    pub rate_ratio: f64,
}

/// Area-contact pressure moments, normalized by the point resultant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointPressureMoments {
    pub peak_pressure_pa: f64,
    pub resultant_n: f64,
    pub second_moment_m2: f64,
}

/// Semiaxes of an elliptic point-contact patch in the principal-curvature
/// frame. The major axis lies along the smaller supplied curvature.
///
/// This is present only for [`NormalPatchLaw::HertzEllipticParaboloid`]. It is
/// not an effective-radius approximation: `patch_radius_m` is an
/// area-equivalent reporting value retained for schema compatibility and must
/// not be used for physical ranking or applicability decisions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EllipticPatchAxes {
    pub semi_major_axis_m: f64,
    pub semi_minor_axis_m: f64,
    pub aspect_ratio: f64,
}

/// Line-contact pressure moments, normalized by the line resultant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinePressureMoments {
    pub peak_pressure_pa: f64,
    pub resultant_n_per_m: f64,
    pub second_moment_m2: f64,
}

/// Point-contact physical compliant normal response.
#[derive(Debug, Clone, PartialEq)]
pub struct PointNormalPatchReceipt {
    pub request_id: ContentHash,
    pub receipt_id: ContentHash,
    pub units: &'static str,
    pub interface_system_id: String,
    pub history_id: String,
    pub input_source_id: String,
    pub authority: InputAuthority,
    pub approach_m: f64,
    pub normal_force_n: f64,
    pub tangent_n_per_m: f64,
    /// Circular-patch radius. For an elliptic receipt this is only the
    /// area-equivalent reporting value `sqrt(a*b)`; it does not enter the
    /// constitutive force, tangent, or applicability admission. Consumers
    /// making a physical comparison must use `elliptic_patch_axes`.
    pub patch_radius_m: f64,
    /// Exact elliptic axes when the admitted law produces an ellipse.
    pub elliptic_patch_axes: Option<EllipticPatchAxes>,
    pub pressure: PointPressureMoments,
    pub reversible_energy_j: f64,
    pub irreversible_work_j: f64,
    pub dissipated_power_w: f64,
    pub ratios: ApplicabilityRatios,
    pub uncertainty: InputUncertainty,
}

/// Line-contact physical response. Every extensive quantity is per unit axial length.
#[derive(Debug, Clone, PartialEq)]
pub struct LineNormalPatchReceipt {
    pub request_id: ContentHash,
    pub receipt_id: ContentHash,
    pub units: &'static str,
    pub interface_system_id: String,
    pub history_id: String,
    pub input_source_id: String,
    pub authority: InputAuthority,
    pub approach_m: f64,
    pub normal_line_load_n_per_m: f64,
    pub tangent_pa: f64,
    pub patch_half_width_m: f64,
    pub pressure: LinePressureMoments,
    pub reversible_energy_j_per_m: f64,
    pub irreversible_work_j_per_m: f64,
    pub dissipated_power_w_per_m: f64,
    pub ratios: ApplicabilityRatios,
    pub uncertainty: InputUncertainty,
}

/// Distinct physical receipts prevent point and line units from being mixed.
#[derive(Debug, Clone, PartialEq)]
pub enum NormalPatchReceipt {
    Point(PointNormalPatchReceipt),
    Line(LineNormalPatchReceipt),
}

impl NormalPatchReceipt {
    /// Canonical identity of the request that produced this response.
    pub const fn request_id(&self) -> ContentHash {
        match self {
            Self::Point(receipt) => receipt.request_id,
            Self::Line(receipt) => receipt.request_id,
        }
    }

    /// Canonical identity of this physical response.
    pub const fn receipt_id(&self) -> ContentHash {
        match self {
            Self::Point(receipt) => receipt.receipt_id,
            Self::Line(receipt) => receipt.receipt_id,
        }
    }
}

/// A nonphysical constraint/barrier marker. It has no force, patch, or energy fields.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConstraintBarrierReceipt {
    pub request_id: ContentHash,
    pub scheme_id: String,
}

/// Total refusal surface for finite-patch law admission/evaluation.
#[derive(Debug, Clone, PartialEq)]
pub enum NormalPatchError {
    MissingIdentity {
        field: &'static str,
    },
    InvalidInput {
        field: &'static str,
    },
    OutsideApplicability {
        ratio: &'static str,
        value: f64,
        limit: f64,
    },
    AdhesionUnsupported {
        adhesion_energy_j_per_m2: f64,
    },
    UnsupportedGeometry {
        geometry: NormalPatchGeometry,
    },
    GeometryLawMismatch {
        geometry: NormalPatchGeometry,
    },
    /// The bounded elliptic-integral solver intentionally refuses aspect
    /// ratios too extreme for its documented numerical envelope.
    EllipticAspectRatioUnsupported {
        curvature_ratio: f64,
        maximum_supported: f64,
    },
    /// A finite in-envelope elliptic-integral solve did not converge.
    EllipticShapeSolveFailure,
    Overflow {
        field: &'static str,
    },
}

impl fmt::Display for NormalPatchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingIdentity { field } => write!(f, "nonblank identity required: {field}"),
            Self::InvalidInput { field } => write!(f, "invalid finite-patch input: {field}"),
            Self::OutsideApplicability {
                ratio,
                value,
                limit,
            } => write!(f, "{ratio}={value} exceeds declared limit {limit}"),
            Self::AdhesionUnsupported {
                adhesion_energy_j_per_m2,
            } => write!(
                f,
                "nonadhesive Hertz/Hunt-Crossley rung refuses adhesion energy {adhesion_energy_j_per_m2}"
            ),
            Self::UnsupportedGeometry { geometry } => write!(
                f,
                "normal-patch law has no toroidal/highly-elliptical rung: {geometry:?}"
            ),
            Self::GeometryLawMismatch { geometry } => write!(
                f,
                "normal-patch law does not match declared local geometry: {geometry:?}"
            ),
            Self::EllipticAspectRatioUnsupported {
                curvature_ratio,
                maximum_supported,
            } => write!(
                f,
                "elliptic Hertz curvature ratio {curvature_ratio} exceeds the documented numerical limit {maximum_supported}"
            ),
            Self::EllipticShapeSolveFailure => write!(
                f,
                "elliptic Hertz aspect-ratio solve did not converge within its declared numerical envelope"
            ),
            Self::Overflow { field } => write!(f, "finite-patch derived value overflowed: {field}"),
        }
    }
}
impl std::error::Error for NormalPatchError {}

impl NormalPatchRequest {
    /// Evaluates one bounded request without mutating solver or history state.
    pub fn evaluate(&self) -> Result<NormalPatchReceipt, NormalPatchError> {
        self.validate()?;
        let request_id = hash_domain(DOMAIN, self.canonical().as_bytes());
        let (radius, modulus, dissipative) = match self.law {
            NormalPatchLaw::HertzSpherePlane {
                effective_radius_m,
                reduced_modulus_pa,
            } => (effective_radius_m, reduced_modulus_pa, None),
            NormalPatchLaw::HertzCylinderPlane {
                effective_radius_m,
                reduced_modulus_pa,
            } => return self.cylinder(request_id, effective_radius_m, reduced_modulus_pa),
            NormalPatchLaw::HertzEllipticParaboloid {
                maximum_principal_curvature_m_inverse,
                minimum_principal_curvature_m_inverse,
                reduced_modulus_pa,
            } => {
                return self.elliptic_paraboloid(
                    request_id,
                    maximum_principal_curvature_m_inverse,
                    minimum_principal_curvature_m_inverse,
                    reduced_modulus_pa,
                    None,
                );
            }
            NormalPatchLaw::HuntCrossleyEllipticParaboloid {
                maximum_principal_curvature_m_inverse,
                minimum_principal_curvature_m_inverse,
                reduced_modulus_pa,
                dissipation_s_per_m,
            } => {
                return self.elliptic_paraboloid(
                    request_id,
                    maximum_principal_curvature_m_inverse,
                    minimum_principal_curvature_m_inverse,
                    reduced_modulus_pa,
                    Some(dissipation_s_per_m),
                );
            }
            NormalPatchLaw::HuntCrossleySphere {
                effective_radius_m,
                reduced_modulus_pa,
                dissipation_s_per_m,
            } => (
                effective_radius_m,
                reduced_modulus_pa,
                Some(dissipation_s_per_m),
            ),
        };
        match self.geometry {
            NormalPatchGeometry::SpherePlane => {
                self.sphere(request_id, radius, modulus, dissipative)
            }
            NormalPatchGeometry::CylinderPlane => Err(NormalPatchError::GeometryLawMismatch {
                geometry: self.geometry,
            }),
            NormalPatchGeometry::EllipticParaboloid => Err(NormalPatchError::GeometryLawMismatch {
                geometry: self.geometry,
            }),
            NormalPatchGeometry::ToroidalOrHighlyElliptical => {
                Err(NormalPatchError::UnsupportedGeometry {
                    geometry: self.geometry,
                })
            }
        }
    }

    fn sphere(
        &self,
        request_id: ContentHash,
        radius: f64,
        modulus: f64,
        dissipative: Option<f64>,
    ) -> Result<NormalPatchReceipt, NormalPatchError> {
        let d = self.indentation_m;
        let a = checked((radius * d).sqrt(), "sphere_patch")?;
        let elastic_force = checked(
            (4.0 / 3.0) * modulus * radius.sqrt() * d.powf(1.5),
            "sphere_force",
        )?;
        let elastic_tangent = if d == 0.0 {
            0.0
        } else {
            checked(2.0 * modulus * (radius * d).sqrt(), "sphere_tangent")?
        };
        let factor = dissipative.map_or(1.0, |c| 1.0 + c * self.indentation_rate_m_per_s);
        if factor < 0.0 {
            return Err(NormalPatchError::OutsideApplicability {
                ratio: "Hunt-Crossley force factor",
                value: factor,
                limit: 0.0,
            });
        }
        let force = checked(elastic_force * factor, "normal_force")?;
        let tangent = checked(elastic_tangent * factor, "normal_tangent")?;
        let p0 = if a == 0.0 {
            0.0
        } else {
            checked(
                1.5 * force / (core::f64::consts::PI * a * a),
                "peak_pressure",
            )?
        };
        let ratios = self.ratios(a, radius, d, p0)?;
        let irreversible_power = dissipative.map_or(0.0, |c| {
            elastic_force * c * self.indentation_rate_m_per_s.powi(2)
        });
        let irreversible_work = checked(irreversible_power * self.step_s, "irreversible_work")?;
        let reversible = checked(0.4 * elastic_force * d, "reversible_energy")?;
        self.point_receipt(
            request_id,
            d,
            force,
            tangent,
            a,
            p0,
            0.4 * a * a,
            reversible,
            irreversible_work,
            irreversible_power,
            ratios,
        )
    }

    fn cylinder(
        &self,
        request_id: ContentHash,
        radius: f64,
        modulus: f64,
    ) -> Result<NormalPatchReceipt, NormalPatchError> {
        let q = self.line_load_n_per_m;
        if q == 0.0 {
            let ratios = self.ratios(0.0, radius, 0.0, 0.0)?;
            return self.line_receipt(
                request_id, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ratios,
            );
        }
        let a = checked(
            (4.0 * q * radius / (core::f64::consts::PI * modulus)).sqrt(),
            "cylinder_patch",
        )?;
        let log = checked(
            (4.0 * radius / a.max(f64::MIN_POSITIVE)).ln(),
            "cylinder_log",
        )?;
        if log <= 0.0 {
            return Err(NormalPatchError::OutsideApplicability {
                ratio: "cylinder logarithmic reach",
                value: log,
                limit: 0.0,
            });
        }
        let approach = checked(
            q * (log + 0.5) / (core::f64::consts::PI * modulus),
            "cylinder_approach",
        )?;
        let tangent = checked(core::f64::consts::PI * modulus / log, "cylinder_tangent")?;
        let p0 = checked(2.0 * q / (core::f64::consts::PI * a), "cylinder_pressure")?;
        let ratios = self.ratios(a, radius, approach, p0)?;
        let energy = checked(
            q * q * (log + 0.75) / (2.0 * core::f64::consts::PI * modulus),
            "cylinder_energy",
        )?;
        self.line_receipt(
            request_id,
            approach,
            q,
            tangent,
            a,
            p0,
            0.25 * a * a,
            energy,
            0.0,
            0.0,
            ratios,
        )
    }

    /// Evaluates the elastic Hertz solution for a two-positive-curvature local
    /// gap. The complete elliptic integrals determine the ellipse aspect ratio
    /// from the supplied curvatures; no scalar effective-radius contact law is
    /// used.
    ///
    /// With `A = k_min / 2`, `B = k_max / 2`, and major/minor semiaxes
    /// `a >= b`, the construction is
    ///
    /// ```text
    /// B/A = (a^2 E(m) - b^2 K(m)) / (b^2 (K(m) - E(m))),
    /// m = 1 - (b/a)^2,
    /// b^2 = d E(m) / ((A + B) K(m)),
    /// F = 2 pi E* b^3 (a/b) (A + B) / (3 E(m)).
    /// ```
    ///
    /// `K` and `E` are complete elliptic integrals. This is the classical
    /// nonadhesive, small-strain half-space solution. It does not establish
    /// that a finite chamfer is locally parabolic, nor does it include
    /// plasticity, adhesion, viscoelasticity, or tangential coupling.
    fn elliptic_paraboloid(
        &self,
        request_id: ContentHash,
        maximum_curvature: f64,
        minimum_curvature: f64,
        modulus: f64,
        dissipative: Option<f64>,
    ) -> Result<NormalPatchReceipt, NormalPatchError> {
        let d = self.indentation_m;
        let shape = elliptic_hertz_shape(maximum_curvature, minimum_curvature)?;
        let curvature_sum = checked(maximum_curvature + minimum_curvature, "curvature_sum")?;
        let b_squared = checked(
            d * shape.complete_e / (0.5 * curvature_sum * shape.complete_k),
            "elliptic_minor_axis_squared",
        )?;
        let b = checked(b_squared.sqrt(), "elliptic_minor_axis")?;
        let a = checked(b * shape.aspect_ratio, "elliptic_major_axis")?;
        let elastic_force = checked(
            core::f64::consts::PI * modulus * b.powi(3) * shape.aspect_ratio * curvature_sum
                / (3.0 * shape.complete_e),
            "elliptic_elastic_force",
        )?;
        let elastic_tangent = if d == 0.0 {
            0.0
        } else {
            checked(1.5 * elastic_force / d, "elliptic_elastic_tangent")?
        };
        let factor = dissipative.map_or(1.0, |c| 1.0 + c * self.indentation_rate_m_per_s);
        if factor < 0.0 {
            return Err(NormalPatchError::OutsideApplicability {
                ratio: "Hunt-Crossley force factor",
                value: factor,
                limit: 0.0,
            });
        }
        let force = checked(elastic_force * factor, "elliptic_normal_force")?;
        let tangent = checked(elastic_tangent * factor, "elliptic_normal_tangent")?;
        let patch_area = checked(core::f64::consts::PI * a * b, "elliptic_patch_area")?;
        let p0 = if patch_area == 0.0 {
            0.0
        } else {
            checked(1.5 * force / patch_area, "elliptic_peak_pressure")?
        };
        // `a` and `1/k_max` are deliberately conservative scales: they make
        // the admitted small-patch and small-strain envelope hold in both
        // principal directions without replacing the constitutive law by an
        // effective-radius approximation.
        let minimum_radius = checked(1.0 / maximum_curvature, "minimum_radius")?;
        let ratios = self.ratios(a, minimum_radius, d, p0)?;
        let axes = EllipticPatchAxes {
            semi_major_axis_m: a,
            semi_minor_axis_m: b,
            aspect_ratio: shape.aspect_ratio,
        };
        let reporting_radius = checked((a * b).sqrt(), "elliptic_reporting_radius")?;
        let reversible = checked(0.4 * elastic_force * d, "elliptic_reversible_energy")?;
        let irreversible_power = dissipative.map_or(0.0, |c| {
            elastic_force * c * self.indentation_rate_m_per_s.powi(2)
        });
        let irreversible_work = checked(
            irreversible_power * self.step_s,
            "elliptic_irreversible_work",
        )?;
        self.elliptic_point_receipt(
            request_id,
            d,
            force,
            tangent,
            reporting_radius,
            axes,
            p0,
            checked((a * a + b * b) / 5.0, "elliptic_second_moment")?,
            reversible,
            irreversible_work,
            irreversible_power,
            ratios,
        )
    }

    fn point_receipt(
        &self,
        request_id: ContentHash,
        approach: f64,
        force: f64,
        tangent: f64,
        patch: f64,
        peak: f64,
        moment: f64,
        reversible: f64,
        irreversible: f64,
        power: f64,
        ratios: ApplicabilityRatios,
    ) -> Result<NormalPatchReceipt, NormalPatchError> {
        self.validate_response(
            approach,
            force,
            tangent,
            patch,
            peak,
            reversible,
            irreversible,
            power,
        )?;
        let authority = self.interface.provenance().authority();
        let receipt_id = self.receipt_id(
            "point",
            request_id,
            approach,
            force,
            tangent,
            patch,
            peak,
            moment,
            reversible,
            irreversible,
            power,
            ratios,
        );
        Ok(NormalPatchReceipt::Point(PointNormalPatchReceipt {
            request_id,
            receipt_id,
            units: POINT_SI_UNITS,
            interface_system_id: self.interface.ordered_system_id().to_owned(),
            history_id: self.interface.history_id().to_owned(),
            input_source_id: self.interface.provenance().source_id().to_owned(),
            authority,
            approach_m: approach,
            normal_force_n: force,
            tangent_n_per_m: tangent,
            patch_radius_m: patch,
            elliptic_patch_axes: None,
            pressure: PointPressureMoments {
                peak_pressure_pa: peak,
                resultant_n: force,
                second_moment_m2: moment,
            },
            reversible_energy_j: reversible,
            irreversible_work_j: irreversible,
            dissipated_power_w: power,
            ratios,
            uncertainty: self.uncertainty,
        }))
    }

    #[allow(clippy::too_many_arguments)]
    fn elliptic_point_receipt(
        &self,
        request_id: ContentHash,
        approach: f64,
        force: f64,
        tangent: f64,
        reporting_radius: f64,
        axes: EllipticPatchAxes,
        peak: f64,
        moment: f64,
        reversible: f64,
        irreversible: f64,
        power: f64,
        ratios: ApplicabilityRatios,
    ) -> Result<NormalPatchReceipt, NormalPatchError> {
        self.validate_response(
            approach,
            force,
            tangent,
            reporting_radius,
            peak,
            reversible,
            irreversible,
            power,
        )?;
        let receipt_id = self.elliptic_receipt_id(
            request_id,
            approach,
            force,
            tangent,
            reporting_radius,
            axes,
            peak,
            moment,
            reversible,
            irreversible,
            power,
            ratios,
        );
        Ok(NormalPatchReceipt::Point(PointNormalPatchReceipt {
            request_id,
            receipt_id,
            units: POINT_SI_UNITS,
            interface_system_id: self.interface.ordered_system_id().to_owned(),
            history_id: self.interface.history_id().to_owned(),
            input_source_id: self.interface.provenance().source_id().to_owned(),
            authority: self.interface.provenance().authority(),
            approach_m: approach,
            normal_force_n: force,
            tangent_n_per_m: tangent,
            patch_radius_m: reporting_radius,
            elliptic_patch_axes: Some(axes),
            pressure: PointPressureMoments {
                peak_pressure_pa: peak,
                resultant_n: force,
                second_moment_m2: moment,
            },
            reversible_energy_j: reversible,
            irreversible_work_j: irreversible,
            dissipated_power_w: power,
            ratios,
            uncertainty: self.uncertainty,
        }))
    }

    fn line_receipt(
        &self,
        request_id: ContentHash,
        approach: f64,
        line_load: f64,
        tangent: f64,
        half_width: f64,
        peak: f64,
        moment: f64,
        reversible: f64,
        irreversible: f64,
        power: f64,
        ratios: ApplicabilityRatios,
    ) -> Result<NormalPatchReceipt, NormalPatchError> {
        self.validate_response(
            approach,
            line_load,
            tangent,
            half_width,
            peak,
            reversible,
            irreversible,
            power,
        )?;
        let receipt_id = self.receipt_id(
            "line",
            request_id,
            approach,
            line_load,
            tangent,
            half_width,
            peak,
            moment,
            reversible,
            irreversible,
            power,
            ratios,
        );
        Ok(NormalPatchReceipt::Line(LineNormalPatchReceipt {
            request_id,
            receipt_id,
            units: LINE_SI_UNITS,
            interface_system_id: self.interface.ordered_system_id().to_owned(),
            history_id: self.interface.history_id().to_owned(),
            input_source_id: self.interface.provenance().source_id().to_owned(),
            authority: self.interface.provenance().authority(),
            approach_m: approach,
            normal_line_load_n_per_m: line_load,
            tangent_pa: tangent,
            patch_half_width_m: half_width,
            pressure: LinePressureMoments {
                peak_pressure_pa: peak,
                resultant_n_per_m: line_load,
                second_moment_m2: moment,
            },
            reversible_energy_j_per_m: reversible,
            irreversible_work_j_per_m: irreversible,
            dissipated_power_w_per_m: power,
            ratios,
            uncertainty: self.uncertainty,
        }))
    }

    fn validate_response(
        &self,
        approach: f64,
        force: f64,
        tangent: f64,
        patch: f64,
        peak: f64,
        reversible: f64,
        irreversible: f64,
        power: f64,
    ) -> Result<(), NormalPatchError> {
        for (value, field) in [
            (approach, "approach"),
            (force, "force"),
            (tangent, "tangent"),
            (patch, "patch"),
            (peak, "pressure"),
            (reversible, "energy"),
            (irreversible, "work"),
            (power, "power"),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(NormalPatchError::Overflow { field });
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn receipt_id(
        &self,
        kind: &str,
        request_id: ContentHash,
        approach: f64,
        force: f64,
        tangent: f64,
        patch: f64,
        peak: f64,
        moment: f64,
        reversible: f64,
        irreversible: f64,
        power: f64,
        ratios: ApplicabilityRatios,
    ) -> ContentHash {
        hash_domain(
            DOMAIN,
            format!(
                concat!(
                    "{}|{}|{}|{:.17e}|{:.17e}|{:.17e}|",
                    "{:.17e}|{:.17e}|{:.17e}|{:.17e}|",
                    "{:.17e}|{:.17e}|{:?}|{:.17e}|{:.17e}|",
                    "{:.17e}|{:?}|{}|{}|{}"
                ),
                kind,
                request_id.to_hex(),
                self.identity.state_id,
                approach,
                force,
                tangent,
                patch,
                peak,
                moment,
                reversible,
                irreversible,
                power,
                ratios,
                self.uncertainty.radius_relative,
                self.uncertainty.modulus_relative,
                self.uncertainty.load_relative,
                self.interface.provenance().authority(),
                self.interface.ordered_system_id(),
                self.interface.history_id(),
                self.interface.provenance().source_id(),
            )
            .as_bytes(),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn elliptic_receipt_id(
        &self,
        request_id: ContentHash,
        approach: f64,
        force: f64,
        tangent: f64,
        reporting_radius: f64,
        axes: EllipticPatchAxes,
        peak: f64,
        moment: f64,
        reversible: f64,
        irreversible: f64,
        power: f64,
        ratios: ApplicabilityRatios,
    ) -> ContentHash {
        hash_domain(
            DOMAIN,
            format!(
                concat!(
                    "elliptic-point|{}|{}|{:.17e}|{:.17e}|{:.17e}|{:.17e}|",
                    "{:.17e}|{:.17e}|{:.17e}|{:.17e}|{:.17e}|{:.17e}|",
                    "{:.17e}|{:.17e}|{:?}|{:.17e}|{:.17e}|{:.17e}|{:?}|{}|{}|{}"
                ),
                request_id.to_hex(),
                self.identity.state_id,
                approach,
                force,
                tangent,
                reporting_radius,
                axes.semi_major_axis_m,
                axes.semi_minor_axis_m,
                axes.aspect_ratio,
                peak,
                moment,
                reversible,
                irreversible,
                power,
                ratios,
                self.uncertainty.radius_relative,
                self.uncertainty.modulus_relative,
                self.uncertainty.load_relative,
                self.interface.provenance().authority(),
                self.interface.ordered_system_id(),
                self.interface.history_id(),
                self.interface.provenance().source_id(),
            )
            .as_bytes(),
        )
    }

    fn ratios(
        &self,
        patch: f64,
        radius: f64,
        approach: f64,
        pressure: f64,
    ) -> Result<ApplicabilityRatios, NormalPatchError> {
        let r = ApplicabilityRatios {
            patch_to_radius: patch / radius,
            strain: approach / radius,
            patch_to_depth: patch / self.applicability.half_space_depth_m,
            patch_to_layer: patch / self.applicability.layer_thickness_m,
            pressure_to_yield: pressure / self.applicability.yield_strength_pa,
            rate_ratio: self.indentation_rate_m_per_s.abs()
                / self.applicability.characteristic_rate_m_per_s,
        };
        for (name, value, limit) in [
            (
                "patch_to_radius",
                r.patch_to_radius,
                self.limits.max_patch_to_radius,
            ),
            ("strain", r.strain, self.limits.max_strain),
            (
                "patch_to_depth",
                r.patch_to_depth,
                self.limits.max_patch_to_depth,
            ),
            (
                "patch_to_layer",
                r.patch_to_layer,
                self.limits.max_patch_to_layer,
            ),
            (
                "pressure_to_yield",
                r.pressure_to_yield,
                self.limits.max_pressure_to_yield,
            ),
            ("rate_ratio", r.rate_ratio, self.limits.max_rate_ratio),
        ] {
            if !value.is_finite() {
                return Err(NormalPatchError::Overflow { field: name });
            }
            if value > limit {
                return Err(NormalPatchError::OutsideApplicability {
                    ratio: name,
                    value,
                    limit,
                });
            }
        }
        Ok(r)
    }

    fn validate(&self) -> Result<(), NormalPatchError> {
        for (s, name) in [
            (&self.identity.model_id, "model_id"),
            (&self.identity.source_id, "source_id"),
            (&self.identity.state_id, "state_id"),
        ] {
            if s.trim().is_empty() {
                return Err(NormalPatchError::MissingIdentity { field: name });
            }
        }
        for (value, name) in [
            (self.indentation_m, "indentation_m"),
            (self.step_s, "step_s"),
            (self.line_load_n_per_m, "line_load_n_per_m"),
            (self.applicability.half_space_depth_m, "half_space_depth_m"),
            (self.applicability.layer_thickness_m, "layer_thickness_m"),
            (self.applicability.yield_strength_pa, "yield_strength_pa"),
            (
                self.applicability.characteristic_rate_m_per_s,
                "characteristic_rate",
            ),
            (self.applicability.temperature_k, "temperature_k"),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(NormalPatchError::InvalidInput { field: name });
            }
        }
        if self.applicability.half_space_depth_m == 0.0
            || self.applicability.layer_thickness_m == 0.0
            || self.applicability.yield_strength_pa == 0.0
            || self.applicability.characteristic_rate_m_per_s == 0.0
        {
            return Err(NormalPatchError::InvalidInput {
                field: "applicability denominator",
            });
        }
        if !self.indentation_rate_m_per_s.is_finite() {
            return Err(NormalPatchError::InvalidInput {
                field: "indentation_rate",
            });
        }
        if self.applicability.adhesion_energy_j_per_m2 != 0.0 {
            return Err(NormalPatchError::AdhesionUnsupported {
                adhesion_energy_j_per_m2: self.applicability.adhesion_energy_j_per_m2,
            });
        }
        if self.applicability.temperature_k < self.limits.min_temperature_k
            || self.applicability.temperature_k > self.limits.max_temperature_k
        {
            return Err(NormalPatchError::OutsideApplicability {
                ratio: "temperature_k",
                value: self.applicability.temperature_k,
                limit: self.limits.max_temperature_k,
            });
        }
        for (value, name) in [
            (self.limits.max_patch_to_radius, "max_patch_to_radius"),
            (self.limits.max_strain, "max_strain"),
            (self.limits.max_patch_to_depth, "max_patch_to_depth"),
            (self.limits.max_patch_to_layer, "max_patch_to_layer"),
            (self.limits.max_pressure_to_yield, "max_pressure_to_yield"),
            (self.limits.max_rate_ratio, "max_rate_ratio"),
            (self.limits.min_temperature_k, "min_temperature_k"),
            (self.limits.max_temperature_k, "max_temperature_k"),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(NormalPatchError::InvalidInput { field: name });
            }
        }
        if self.limits.min_temperature_k > self.limits.max_temperature_k {
            return Err(NormalPatchError::InvalidInput {
                field: "temperature limits",
            });
        }
        for (v, n) in [
            (self.uncertainty.radius_relative, "radius_uncertainty"),
            (self.uncertainty.modulus_relative, "modulus_uncertainty"),
            (self.uncertainty.load_relative, "load_uncertainty"),
        ] {
            if !v.is_finite() || v < 0.0 {
                return Err(NormalPatchError::InvalidInput { field: n });
            }
        }
        match self.law {
            NormalPatchLaw::HertzSpherePlane {
                effective_radius_m,
                reduced_modulus_pa,
            }
            | NormalPatchLaw::HuntCrossleySphere {
                effective_radius_m,
                reduced_modulus_pa,
                ..
            }
            | NormalPatchLaw::HertzCylinderPlane {
                effective_radius_m,
                reduced_modulus_pa,
            } => {
                if !(effective_radius_m.is_finite()
                    && effective_radius_m > 0.0
                    && reduced_modulus_pa.is_finite()
                    && reduced_modulus_pa > 0.0)
                {
                    return Err(NormalPatchError::InvalidInput {
                        field: "curvature_or_modulus",
                    });
                }
            }
            NormalPatchLaw::HertzEllipticParaboloid {
                maximum_principal_curvature_m_inverse,
                minimum_principal_curvature_m_inverse,
                reduced_modulus_pa,
            }
            | NormalPatchLaw::HuntCrossleyEllipticParaboloid {
                maximum_principal_curvature_m_inverse,
                minimum_principal_curvature_m_inverse,
                reduced_modulus_pa,
                ..
            } => {
                if !(maximum_principal_curvature_m_inverse.is_finite()
                    && maximum_principal_curvature_m_inverse > 0.0
                    && minimum_principal_curvature_m_inverse.is_finite()
                    && minimum_principal_curvature_m_inverse > 0.0
                    && reduced_modulus_pa.is_finite()
                    && reduced_modulus_pa > 0.0)
                {
                    return Err(NormalPatchError::InvalidInput {
                        field: "principal_curvature_or_modulus",
                    });
                }
                if minimum_principal_curvature_m_inverse > maximum_principal_curvature_m_inverse {
                    return Err(NormalPatchError::InvalidInput {
                        field: "principal_curvature_order",
                    });
                }
            }
        }
        if let NormalPatchLaw::HuntCrossleySphere {
            dissipation_s_per_m,
            ..
        }
        | NormalPatchLaw::HuntCrossleyEllipticParaboloid {
            dissipation_s_per_m,
            ..
        } = self.law
        {
            if !dissipation_s_per_m.is_finite() || dissipation_s_per_m < 0.0 {
                return Err(NormalPatchError::InvalidInput {
                    field: "dissipation",
                });
            }
        }
        match (self.law, self.geometry) {
            (_, NormalPatchGeometry::ToroidalOrHighlyElliptical) => {
                return Err(NormalPatchError::UnsupportedGeometry {
                    geometry: self.geometry,
                });
            }
            (NormalPatchLaw::HertzCylinderPlane { .. }, NormalPatchGeometry::CylinderPlane)
            | (
                NormalPatchLaw::HertzSpherePlane { .. } | NormalPatchLaw::HuntCrossleySphere { .. },
                NormalPatchGeometry::SpherePlane,
            )
            | (
                NormalPatchLaw::HertzEllipticParaboloid { .. },
                NormalPatchGeometry::EllipticParaboloid,
            )
            | (
                NormalPatchLaw::HuntCrossleyEllipticParaboloid { .. },
                NormalPatchGeometry::EllipticParaboloid,
            ) => {}
            _ => {
                return Err(NormalPatchError::GeometryLawMismatch {
                    geometry: self.geometry,
                });
            }
        }
        Ok(())
    }

    fn canonical(&self) -> String {
        format!(
            concat!(
                "v1|{}|{}|{}|{:?}|{:?}|{:.17e}|{:.17e}|{:.17e}|{:.17e}|",
                "{:.17e}|{:.17e}|{:.17e}|{:.17e}|{:.17e}|{:.17e}|",
                "{:.17e}|{:.17e}|{:.17e}|{:.17e}|{:.17e}|{:.17e}|",
                "{:.17e}|{:.17e}|{:.17e}|{:.17e}|{:.17e}|{}|{}|{}|{:?}"
            ),
            self.identity.model_id,
            self.identity.source_id,
            self.identity.state_id,
            self.law,
            self.geometry,
            self.indentation_m,
            self.indentation_rate_m_per_s,
            self.step_s,
            self.line_load_n_per_m,
            self.applicability.half_space_depth_m,
            self.applicability.layer_thickness_m,
            self.applicability.yield_strength_pa,
            self.applicability.characteristic_rate_m_per_s,
            self.applicability.temperature_k,
            self.applicability.adhesion_energy_j_per_m2,
            self.limits.max_patch_to_radius,
            self.limits.max_strain,
            self.limits.max_patch_to_depth,
            self.limits.max_patch_to_layer,
            self.limits.max_pressure_to_yield,
            self.limits.max_rate_ratio,
            self.limits.min_temperature_k,
            self.limits.max_temperature_k,
            self.uncertainty.radius_relative,
            self.uncertainty.modulus_relative,
            self.uncertainty.load_relative,
            self.interface.ordered_system_id(),
            self.interface.history_id(),
            self.interface.provenance().source_id(),
            self.interface.provenance().authority(),
        )
    }
}

fn checked(value: f64, field: &'static str) -> Result<f64, NormalPatchError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(NormalPatchError::Overflow { field })
    }
}

/// Numerical shape factors for the classical elliptic Hertz solution.
#[derive(Debug, Clone, Copy)]
struct EllipticHertzShape {
    aspect_ratio: f64,
    complete_k: f64,
    complete_e: f64,
}

/// Smallest permitted squared minor-to-major axis ratio. The associated
/// curvature-ratio threshold is derived from the elliptic equation rather
/// than guessed from an effective radius. The limit keeps `1 - m` separated
/// from the logarithmic singular endpoint of the complete integral in f64.
const MIN_MINOR_TO_MAJOR_SQUARED: f64 = 1.0e-12;
const ELLIPTIC_BISECTION_STEPS: usize = 80;

fn elliptic_hertz_shape(
    maximum_curvature: f64,
    minimum_curvature: f64,
) -> Result<EllipticHertzShape, NormalPatchError> {
    if maximum_curvature.to_bits() == minimum_curvature.to_bits() {
        return Ok(EllipticHertzShape {
            aspect_ratio: 1.0,
            complete_k: core::f64::consts::FRAC_PI_2,
            complete_e: core::f64::consts::FRAC_PI_2,
        });
    }
    let target = minimum_curvature / maximum_curvature;
    let (low_k, low_e) = complete_elliptic_integrals(1.0 - MIN_MINOR_TO_MAJOR_SQUARED)?;
    let minimum_target = elliptic_curvature_ratio(MIN_MINOR_TO_MAJOR_SQUARED, low_k, low_e)?;
    if target < minimum_target {
        return Err(NormalPatchError::EllipticAspectRatioUnsupported {
            curvature_ratio: maximum_curvature / minimum_curvature,
            maximum_supported: 1.0 / minimum_target,
        });
    }

    let mut lower = MIN_MINOR_TO_MAJOR_SQUARED;
    let mut upper = 1.0;
    for _ in 0..ELLIPTIC_BISECTION_STEPS {
        let middle = 0.5 * (lower + upper);
        let (k, e) = complete_elliptic_integrals(1.0 - middle)?;
        let recovered = elliptic_curvature_ratio(middle, k, e)?;
        if recovered < target {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    let minor_to_major_squared = 0.5 * (lower + upper);
    let (complete_k, complete_e) = complete_elliptic_integrals(1.0 - minor_to_major_squared)?;
    let recovered = elliptic_curvature_ratio(minor_to_major_squared, complete_k, complete_e)?;
    if (recovered - target).abs() > 1.0e-12 * target.max(MIN_MINOR_TO_MAJOR_SQUARED) {
        return Err(NormalPatchError::EllipticShapeSolveFailure);
    }
    Ok(EllipticHertzShape {
        aspect_ratio: 1.0 / minor_to_major_squared.sqrt(),
        complete_k,
        complete_e,
    })
}

/// Returns `k_min / k_max` from `q = (b/a)^2` and complete elliptic
/// integrals at `m = 1 - q`.
fn elliptic_curvature_ratio(
    minor_to_major_squared: f64,
    complete_k: f64,
    complete_e: f64,
) -> Result<f64, NormalPatchError> {
    if minor_to_major_squared >= 1.0 {
        return Ok(1.0);
    }
    checked(
        (complete_k - complete_e) / (complete_e / minor_to_major_squared - complete_k),
        "elliptic_curvature_ratio",
    )
}

/// Complete elliptic integrals `K(m)` and `E(m)` for `0 <= m < 1`, evaluated
/// by the arithmetic-geometric mean. The series for `E` is the companion AGM
/// identity, so both values share a deterministic, cancellation-free kernel.
fn complete_elliptic_integrals(parameter_m: f64) -> Result<(f64, f64), NormalPatchError> {
    if !(parameter_m.is_finite() && (0.0..1.0).contains(&parameter_m)) {
        return Err(NormalPatchError::EllipticShapeSolveFailure);
    }
    if parameter_m == 0.0 {
        return Ok((core::f64::consts::FRAC_PI_2, core::f64::consts::FRAC_PI_2));
    }
    let mut arithmetic = 1.0;
    let mut geometric = (1.0 - parameter_m).sqrt();
    let mut sum = 0.5 * parameter_m;
    let mut weight = 1.0;
    for _ in 0..32 {
        let difference = 0.5 * (arithmetic - geometric);
        let next_arithmetic = 0.5 * (arithmetic + geometric);
        let next_geometric = (arithmetic * geometric).sqrt();
        sum = checked(sum + weight * difference * difference, "elliptic_e_series")?;
        arithmetic = next_arithmetic;
        geometric = next_geometric;
        if difference.abs() <= 8.0 * f64::EPSILON * arithmetic {
            let complete_k = checked(core::f64::consts::FRAC_PI_2 / arithmetic, "elliptic_k")?;
            let complete_e = checked(complete_k * (1.0 - sum), "elliptic_e")?;
            return Ok((complete_k, complete_e));
        }
        weight *= 2.0;
    }
    Err(NormalPatchError::EllipticShapeSolveFailure)
}

#[cfg(test)]
mod tests {
    use super::{complete_elliptic_integrals, elliptic_hertz_shape};

    #[test]
    fn elliptic_agm_matches_reference_values_at_half_parameter() {
        let (k, e) = complete_elliptic_integrals(0.5).unwrap();
        assert!((k - 1.854_074_677_301_371_9).abs() < 2.0e-15);
        assert!((e - 1.350_643_881_047_675_5).abs() < 2.0e-15);
    }

    #[test]
    fn elliptic_shape_inverts_the_principal_curvature_ratio() {
        let shape = elliptic_hertz_shape(100.0, 25.0).unwrap();
        let q = 1.0 / shape.aspect_ratio.powi(2);
        let recovered =
            (shape.complete_k - shape.complete_e) / (shape.complete_e / q - shape.complete_k);
        assert!((recovered - 0.25).abs() < 3.0e-13);
    }
}

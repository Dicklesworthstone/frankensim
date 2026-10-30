//! First-order deformation loss of a viscoelastic sphere on a rigid plane.
//!
//! The corrected surface-stress relation in Brilliantov, Pimenova and Goldobin,
//! <https://arxiv.org/abs/1410.3289>, equations (54)–(56), is
//! `p_1 = A d_t p_0`, with standard shear/bulk viscosities `eta_s`, `eta_b`.
//! Its positive coefficient equals the low-frequency logarithmic derivative
//! of the isotropic reduced modulus `E* = 4 G (3 K + G) / (3 K + 4 G)`.
//!
//! For an advected Hertz pressure, `d_t p_0 = -V d_x p_0`. Its first moment,
//! integrated by parts, gives `M = -A N R omega`. This rolling inference uses
//! the same pressure-moment argument as Brilliantov and Pöschel (1998),
//! <https://doi.org/10.1209/epl/i1998-00281-7>, but NOT that paper's older
//! nonlinear-viscosity coefficient. The corrected coefficient is linear and
//! monotone in both passive viscosities and remains regular at `nu = 0`.
//!
//! This is a first-order, small-contact, normal-pressure-only Kelvin–Voigt
//! approximation with negligible inertia. It does not resolve the displaced
//! contact boundary, higher-order fields, the traction needed to enforce
//! rolling, or arbitrary-speed viscoelastic response. The original Hertz
//! owner checks the supplied normal patch. No fitted rolling coefficient,
//! independent normal-contact equation, sliding force or spin loss is added.

use super::{
    PatchCurvature, RollingKinematics, RollingLossApplicability, RollingLossError,
    RollingPatchReceipt, finite, nonblank, nonnegative, positive, validate_interface,
};
use crate::{HertzSpherePlane, InputAuthority, InterfaceSystemRef};

pub mod path;

/// Published constitutive form; the material source remains separately named.
pub const VISCOELASTIC_SPHERE_LAW_ID: &str = "first-order Kelvin-Voigt sphere rolling v1";

/// The sphere carrying the viscoelastic bulk loss in the ordered interface.
/// The other surface is the rigid plane for this constitutive approximation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeformingSphereSurface {
    /// The first ordered surface is the deforming sphere.
    SurfaceA,
    /// The second ordered surface is the deforming sphere.
    SurfaceB,
}

/// Explicit sphere geometry and isotropic Kelvin–Voigt material constants.
///
/// The convention is `sigma = 2 G dev(epsilon) + K tr(epsilon) I`
/// `+ 2 eta_s dev(epsilon_dot) + eta_b tr(epsilon_dot) I`.
/// Thus `eta_s` is shear viscosity and `eta_b` is bulk viscosity, both in Pa s.
/// These are material parameters, not a rolling coefficient fitted to a body
/// trajectory. Constant values retain their authored temperature/frequency box.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViscoelasticSphereParameters {
    /// Which ordered body deforms and generates the bulk heat.
    pub sphere_surface: DeformingSphereSurface,
    /// Actual sphere radius [m], not an inferred equivalent radius.
    pub radius_m: f64,
    /// Young modulus of the sphere [Pa].
    pub young_modulus_pa: f64,
    /// Poisson ratio, strictly between -1 and 1/2.
    pub poisson_ratio: f64,
    /// Sphere bulk density [kg/m³], used for the quasistatic wave-speed check.
    pub density_kg_m3: f64,
    /// Standard shear viscosity `eta_s` [Pa s]; corrected paper `eta_1`.
    pub shear_viscosity_pa_s: f64,
    /// Bulk viscosity `eta_b` [Pa s]; corrected paper `eta_2`.
    pub bulk_viscosity_pa_s: f64,
}

/// Explicit caller ceilings on the small parameters of the asymptotic model.
/// Every ceiling lies strictly between zero and one. They are applicability
/// choices, not independently certified bounds on omitted higher-order terms.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QuasistaticSphereLimits {
    /// Maximum Hertz indentation divided by sphere radius.
    pub maximum_indentation_radius_ratio: f64,
    /// Maximum rolling speed divided by transverse elastic wave speed.
    pub maximum_speed_wave_ratio: f64,
    /// Maximum material relaxation time divided by `indentation / |speed|`.
    /// This is the stricter passage time used in the original derivation.
    pub maximum_relaxation_ratio: f64,
}

impl QuasistaticSphereLimits {
    fn validate(self) -> Result<(), RollingLossError> {
        for (field, value) in [
            (
                "maximum_indentation_radius_ratio",
                self.maximum_indentation_radius_ratio,
            ),
            ("maximum_speed_wave_ratio", self.maximum_speed_wave_ratio),
            ("maximum_relaxation_ratio", self.maximum_relaxation_ratio),
        ] {
            if !value.is_finite() || value <= 0.0 || value >= 1.0 {
                return Err(RollingLossError::InvalidInput { field });
            }
        }
        Ok(())
    }
}

/// Immutable source parameters and their material-derived rolling time [s].
#[derive(Debug, Clone, PartialEq)]
pub struct ViscoelasticSphereRollingCard {
    model_id: String,
    source_id: String,
    authority: InputAuthority,
    parameters: ViscoelasticSphereParameters,
    limits: QuasistaticSphereLimits,
    applicability: RollingLossApplicability,
    rolling_time_s: f64,
    maximum_material_time_s: f64,
    shear_wave_speed_m_s: f64,
    reduced_modulus_pa: f64,
}

/// Physical rolling response and the actual small parameters admitted for it.
/// All heat is generated in the card's deforming sphere. Normal Hertz storage
/// belongs to the normal contact owner; this steady rolling law adds no store.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViscoelasticSphereResponse {
    /// Torque opposing the declared rolling angular velocity [N m].
    pub rolling_moment_nm: f64,
    /// Nonnegative bulk viscous power [W], before heat transport.
    pub dissipated_power_w: f64,
    /// Original admitted Hertz indentation [m].
    pub indentation_m: f64,
    /// Indentation/radius; zero for a released patch.
    pub indentation_radius_ratio: f64,
    /// Rolling speed/transverse wave speed; zero for a released patch.
    pub speed_wave_ratio: f64,
    /// Maximum material relaxation time times speed/indentation.
    /// Zero for a released patch, which has no active passage-time condition.
    pub relaxation_ratio: f64,
}

impl ViscoelasticSphereRollingCard {
    /// Admit the material-derived first-order sphere law.
    ///
    /// Shear and bulk viscosities are independently nonnegative. The corrected
    /// perturbative coefficient remains valid for auxetic elastic materials,
    /// zero Poisson ratio, and either zero viscosity. Both zero is exactly
    /// elastic. No singular rescaling or inferred material coefficient is used.
    pub fn new(
        model_id: impl Into<String>,
        source_id: impl Into<String>,
        authority: InputAuthority,
        parameters: ViscoelasticSphereParameters,
        limits: QuasistaticSphereLimits,
        applicability: RollingLossApplicability,
    ) -> Result<Self, RollingLossError> {
        let model_id = model_id.into();
        let source_id = source_id.into();
        nonblank(&model_id, "viscoelastic_sphere_model_id")?;
        nonblank(&source_id, "viscoelastic_sphere_source_id")?;
        limits.validate()?;
        let p = parameters;
        positive(p.radius_m, "sphere_radius_m")?;
        positive(p.young_modulus_pa, "sphere_young_modulus_pa")?;
        positive(p.density_kg_m3, "sphere_density_kg_m3")?;
        if !p.poisson_ratio.is_finite() || p.poisson_ratio <= -1.0 || p.poisson_ratio >= 0.5 {
            return Err(invalid("sphere_poisson_ratio"));
        }
        nonnegative(p.shear_viscosity_pa_s, "sphere_shear_viscosity_pa_s")?;
        nonnegative(p.bulk_viscosity_pa_s, "sphere_bulk_viscosity_pa_s")?;
        let nu = p.poisson_ratio;
        let shear = p.young_modulus_pa / (2.0 * (1.0 + nu));
        let bulk = p.young_modulus_pa / (3.0 * (1.0 - 2.0 * nu));
        let reduced_modulus_pa = p.young_modulus_pa / ((1.0 - nu) * (1.0 + nu));
        positive(shear, "sphere_shear_modulus_pa")?;
        positive(bulk, "sphere_bulk_modulus_pa")?;
        positive(reduced_modulus_pa, "sphere_reduced_modulus_pa")?;
        let shear_wave_speed_m_s = (shear / p.density_kg_m3).sqrt();
        positive(shear_wave_speed_m_s, "sphere_shear_wave_speed_m_s")?;
        let maximum_material_time_s =
            (p.shear_viscosity_pa_s / shear).max(p.bulk_viscosity_pa_s / bulk);
        nonnegative(maximum_material_time_s, "sphere_material_relaxation_time_s")?;
        // A = (1+nu)/(E(1-nu)) *
        //     [4 eta_s (1-nu+nu^2)/3 + eta_b (1-2nu)^2].
        // Equivalently A = w_s eta_s/G + w_b eta_b/K, where both
        // weights are positive and sum to one. This form avoids cancellation
        // and preserves proportional viscosity A=tau at nu=0 and for auxetics.
        let bulk_weight = ((1.0 + nu) * (1.0 - 2.0 * nu)) / (3.0 * (1.0 - nu));
        let shear_weight = (2.0 * (1.0 - nu + nu * nu)) / (3.0 * (1.0 - nu));
        let rolling_time_s = shear_weight * (p.shear_viscosity_pa_s / shear)
            + bulk_weight * (p.bulk_viscosity_pa_s / bulk);
        nonnegative(rolling_time_s, "derived_sphere_rolling_time_s")?;
        if (p.shear_viscosity_pa_s > 0.0 || p.bulk_viscosity_pa_s > 0.0) && rolling_time_s == 0.0 {
            return Err(invalid("unrepresentable_sphere_rolling_time_s"));
        }
        Ok(Self {
            model_id,
            source_id,
            authority,
            parameters,
            limits,
            applicability,
            rolling_time_s,
            maximum_material_time_s,
            shear_wave_speed_m_s,
            reduced_modulus_pa,
        })
    }

    /// Authored physical-model identity.
    #[must_use]
    pub fn model_id(&self) -> &str {
        &self.model_id
    }
    /// Original material source identity.
    #[must_use]
    pub fn source_id(&self) -> &str {
        &self.source_id
    }
    /// Original source authority ceiling.
    #[must_use]
    pub const fn authority(&self) -> InputAuthority {
        self.authority
    }
    /// Original sphere geometry and material constants.
    #[must_use]
    pub const fn parameters(&self) -> ViscoelasticSphereParameters {
        self.parameters
    }
    /// Explicit approximation limits.
    #[must_use]
    pub const fn limits(&self) -> QuasistaticSphereLimits {
        self.limits
    }
    /// Original constant-coefficient temperature/frequency domain.
    #[must_use]
    pub const fn applicability(&self) -> RollingLossApplicability {
        self.applicability
    }
    /// Material-derived coefficient in `M = -k_roll N R omega` [s].
    #[must_use]
    pub const fn rolling_time_s(&self) -> f64 {
        self.rolling_time_s
    }

    pub(super) fn validate(&self) -> Result<(), RollingLossError> {
        // All fields are private and construction admits the full material.
        nonblank(&self.model_id, "viscoelastic_sphere_model_id")?;
        nonblank(&self.source_id, "viscoelastic_sphere_source_id")?;
        self.limits.validate()
    }

    /// Evaluate a held, steady rolling component on the original Hertz patch.
    ///
    /// Both principal relative curvatures must equal the actual sphere's 1/R;
    /// an equivalent-radius approximation refuses. Area and load must agree
    /// with the shared Hertz owner for this sphere and a rigid plane. Signed
    /// contour speed must equal `R * rolling_rate` within arithmetic roundoff.
    /// Spin is retained but contributes no torque or heat here. Surface sliding,
    /// normal evolution and body dynamics remain owned by separate consumers.
    pub fn response(
        &self,
        patch: &RollingPatchReceipt,
        interface: &InterfaceSystemRef,
        kinematics: RollingKinematics,
    ) -> Result<ViscoelasticSphereResponse, RollingLossError> {
        self.validate()?;
        patch.validate()?;
        validate_interface(interface)?;
        self.applicability.require_contains(kinematics)?;
        let radius = self.parameters.radius_m;
        let curvature = 1.0 / radius;
        positive(curvature, "sphere_curvature_per_m")?;
        let PatchCurvature::Principal {
            first_per_m,
            second_per_m,
        } = patch.curvature
        else {
            return Err(invalid(
                "sphere_rolling_requires_actual_principal_curvatures",
            ));
        };
        same(first_per_m, curvature, "sphere_first_principal_curvature")?;
        same(second_per_m, curvature, "sphere_second_principal_curvature")?;
        let speed = radius * kinematics.rolling_rate_rad_s;
        same(
            kinematics.contour_speed_mps,
            speed,
            "sphere_no_slip_speed_mapping",
        )?;
        let a2 = patch.contact_area_m2 / core::f64::consts::PI;
        let indentation_m = a2 / radius;
        nonnegative(indentation_m, "sphere_hertz_indentation_m")?;
        let reference = HertzSpherePlane {
            effective_radius: radius,
            reduced_modulus: self.reduced_modulus_pa,
        }
        .response(interface, indentation_m)
        .map_err(|_| invalid("sphere_rolling_original_hertz_patch"))?;
        same(
            patch.normal_load_n,
            reference.normal_force_n,
            "sphere_rolling_hertz_load_and_area",
        )?;
        if patch.normal_load_n == 0.0 {
            return Ok(ViscoelasticSphereResponse {
                rolling_moment_nm: 0.0,
                dissipated_power_w: 0.0,
                indentation_m: 0.0,
                indentation_radius_ratio: 0.0,
                speed_wave_ratio: 0.0,
                relaxation_ratio: 0.0,
            });
        }
        positive(indentation_m, "loaded_sphere_indentation_m")?;
        let indentation_radius_ratio = indentation_m / radius;
        let speed_wave_ratio = speed.abs() / self.shear_wave_speed_m_s;
        let relaxation_ratio = self.maximum_material_time_s * (speed.abs() / indentation_m);
        require_limit(
            indentation_radius_ratio,
            self.limits.maximum_indentation_radius_ratio,
            "sphere_indentation_radius_ratio",
        )?;
        require_limit(
            speed_wave_ratio,
            self.limits.maximum_speed_wave_ratio,
            "sphere_speed_wave_ratio",
        )?;
        require_limit(
            relaxation_ratio,
            self.limits.maximum_relaxation_ratio,
            "sphere_relaxation_ratio",
        )?;
        let rolling_moment_nm = -((self.rolling_time_s * speed) * patch.normal_load_n);
        let dissipated_power_w = -rolling_moment_nm * kinematics.rolling_rate_rad_s;
        finite(rolling_moment_nm, "viscoelastic_sphere_rolling_moment_nm")?;
        nonnegative(dissipated_power_w, "viscoelastic_sphere_dissipated_power_w")?;
        Ok(ViscoelasticSphereResponse {
            rolling_moment_nm,
            dissipated_power_w,
            indentation_m,
            indentation_radius_ratio,
            speed_wave_ratio,
            relaxation_ratio,
        })
    }
}

fn invalid(field: &'static str) -> RollingLossError {
    RollingLossError::InvalidInput { field }
}

fn same(actual: f64, expected: f64, field: &'static str) -> Result<(), RollingLossError> {
    finite(actual, field)?;
    finite(expected, field)?;
    let scale = actual.abs().max(expected.abs()).max(f64::MIN_POSITIVE);
    if (actual - expected).abs() > 512.0 * f64::EPSILON * scale {
        return Err(invalid(field));
    }
    Ok(())
}

fn require_limit(value: f64, maximum: f64, field: &'static str) -> Result<(), RollingLossError> {
    nonnegative(value, field)?;
    if value > maximum {
        return Err(RollingLossError::OutsideApplicability {
            field,
            value,
            minimum: 0.0,
            maximum,
        });
    }
    Ok(())
}

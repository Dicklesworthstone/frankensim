//! Validated spectral dielectric optics for the frontier path tracer.
//!
//! The types in this module make the physical conventions explicit:
//! wavelengths are vacuum nanometres, Cauchy coefficients use micrometres,
//! extinction is in inverse metres, directions point away from the surface,
//! and sampled transmission uses radiance transport.  Presets are explicitly
//! representative rather than measured assets.

use fs_geom::Vec3;
use fs_math::det;

use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN, LiftedSpectrum, lift_rgb};

const PI: f64 = core::f64::consts::PI;
const MIN_ROUGHNESS_ALPHA: f64 = 1.0e-4;
const MAX_ROUGHNESS_ALPHA: f64 = 1.0;
const MAX_GLASS_IOR: f64 = 3.0;
const UNIT_TOLERANCE: f64 = 2.0e-10;

/// Refusal from validated dielectric construction or analytic evaluation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DielectricError {
    /// An index-of-refraction parameter was non-finite or outside the admitted
    /// positive glass domain.
    InvalidIor,
    /// A reference transmittance, extinction coefficient, or distance was
    /// non-finite or outside the Beer-Lambert domain.
    InvalidAbsorption,
    /// A roughness value was outside the admitted GGX interval.
    InvalidRoughness,
    /// A wavelength was outside the tracer's declared visible interval.
    InvalidWavelength,
    /// A direction or normal was non-finite, non-unit, or oriented outside the
    /// interface convention.
    InvalidDirection,
    /// A scalar interface parameter was non-finite or an absolute index was
    /// outside the admitted `[1, 3]` dielectric domain.
    InvalidInterface,
}

impl core::fmt::Display for DielectricError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter.write_str(match self {
            Self::InvalidIor => "invalid dielectric index of refraction",
            Self::InvalidAbsorption => "invalid Beer-Lambert absorption",
            Self::InvalidRoughness => "invalid dielectric GGX roughness",
            Self::InvalidWavelength => "wavelength is outside the visible tracer interval",
            Self::InvalidDirection => "invalid dielectric direction frame",
            Self::InvalidInterface => "invalid dielectric interface parameters",
        })
    }
}

impl core::error::Error for DielectricError {}

/// A positive Cauchy dispersion law
/// `n(lambda_um) = A + B/lambda_um^2 + C/lambda_um^4`.
///
/// The admitted coefficients are nonnegative and the resulting visible-band
/// index must stay in `[1, 3]`.  That bounded domain covers ordinary optical
/// dielectrics used by this product scene while rejecting metals and malformed
/// fits that require a different optical model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CauchyIor {
    a: f64,
    b_um2: f64,
    c_um4: f64,
}

impl CauchyIor {
    /// Construct and validate Cauchy coefficients expressed in micrometres.
    pub fn try_new(a: f64, b_um2: f64, c_um4: f64) -> Result<Self, DielectricError> {
        if !a.is_finite()
            || !b_um2.is_finite()
            || !c_um4.is_finite()
            || a < 1.0
            || b_um2 < 0.0
            || c_um4 < 0.0
        {
            return Err(DielectricError::InvalidIor);
        }
        let candidate = Self {
            a,
            b_um2: if b_um2 == 0.0 { 0.0 } else { b_um2 },
            c_um4: if c_um4 == 0.0 { 0.0 } else { c_um4 },
        };
        let blue = candidate.eval_unchecked(LAMBDA_MIN);
        let red = candidate.eval_unchecked(LAMBDA_MAX);
        if !blue.is_finite() || !red.is_finite() || red < 1.0 || blue > MAX_GLASS_IOR {
            return Err(DielectricError::InvalidIor);
        }
        Ok(candidate)
    }

    /// Construct a nondispersive visible-band index.
    pub fn try_constant(index: f64) -> Result<Self, DielectricError> {
        Self::try_new(index, 0.0, 0.0)
    }

    /// Evaluate the phase index at a visible vacuum wavelength in nanometres.
    pub fn eval(self, wavelength_nm: f64) -> Result<f64, DielectricError> {
        validate_wavelength(wavelength_nm)?;
        Ok(self.eval_unchecked(wavelength_nm))
    }

    /// Canonical `(A, B_um2, C_um4)` parameters for content identity.
    #[must_use]
    pub const fn coefficients(self) -> [f64; 3] {
        [self.a, self.b_um2, self.c_um4]
    }

    /// Whether refraction direction varies over the visible interval.
    #[must_use]
    pub fn is_dispersive(self) -> bool {
        self.b_um2.to_bits() != 0.0_f64.to_bits() || self.c_um4.to_bits() != 0.0_f64.to_bits()
    }

    pub(crate) fn eval_unchecked(self, wavelength_nm: f64) -> f64 {
        let wavelength_um = wavelength_nm * 1.0e-3;
        let inverse_square = 1.0 / (wavelength_um * wavelength_um);
        self.a + self.b_um2 * inverse_square + self.c_um4 * inverse_square * inverse_square
    }
}

/// Canonical source parameters retained by Beer-Lambert absorption.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BeerLambertParameters {
    /// Zero extinction at every wavelength.
    Clear,
    /// Constant spectral extinction in inverse metres.
    Constant {
        /// Nonnegative extinction coefficient (1/m).
        extinction_per_m: f64,
    },
    /// A bounded lifted spectrum defining transmittance at one path length.
    ReferenceRgb {
        /// Declared linear-RGB reference transmittance.
        linear_rgb: [f64; 3],
        /// Positive reference path length (m).
        distance_m: f64,
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum AbsorptionModel {
    Clear,
    Constant(f64),
    ReferenceRgb {
        spectrum: LiftedSpectrum,
        linear_rgb: [f64; 3],
        distance_m: f64,
    },
}

/// Validated homogeneous spectral absorption under Beer-Lambert attenuation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BeerLambertAbsorption {
    model: AbsorptionModel,
}

impl BeerLambertAbsorption {
    /// A medium with no absorption.
    pub const CLEAR: Self = Self {
        model: AbsorptionModel::Clear,
    };

    /// Construct wavelength-independent extinction in inverse metres.
    pub fn try_constant(extinction_per_m: f64) -> Result<Self, DielectricError> {
        if !extinction_per_m.is_finite() || extinction_per_m < 0.0 {
            return Err(DielectricError::InvalidAbsorption);
        }
        if extinction_per_m == 0.0 {
            return Ok(Self::CLEAR);
        }
        Ok(Self {
            model: AbsorptionModel::Constant(extinction_per_m),
        })
    }

    /// Construct a smooth positive extinction spectrum from the declared
    /// linear-RGB transmittance at `distance_m`.
    ///
    /// Exact `[1, 1, 1]` selects [`Self::CLEAR`].  Other components must lie
    /// in `(0, 1]`; clamping is intentionally forbidden.
    pub fn try_from_rgb_transmittance(
        linear_rgb: [f64; 3],
        distance_m: f64,
    ) -> Result<Self, DielectricError> {
        if !distance_m.is_finite()
            || distance_m <= 0.0
            || linear_rgb
                .iter()
                .any(|value| !value.is_finite() || *value <= 0.0 || *value > 1.0)
        {
            return Err(DielectricError::InvalidAbsorption);
        }
        if linear_rgb
            .iter()
            .all(|value| value.to_bits() == 1.0_f64.to_bits())
        {
            return Ok(Self::CLEAR);
        }
        Ok(Self {
            model: AbsorptionModel::ReferenceRgb {
                spectrum: lift_rgb(linear_rgb),
                linear_rgb,
                distance_m,
            },
        })
    }

    /// Spectral extinction coefficient in inverse metres.
    pub fn extinction_per_m(self, wavelength_nm: f64) -> Result<f64, DielectricError> {
        validate_wavelength(wavelength_nm)?;
        Ok(match self.model {
            AbsorptionModel::Clear => 0.0,
            AbsorptionModel::Constant(extinction) => extinction,
            AbsorptionModel::ReferenceRgb {
                spectrum,
                distance_m,
                ..
            } => -det::ln(spectrum.eval(wavelength_nm)) / distance_m,
        })
    }

    /// Transmittance after a nonnegative path length in metres.
    pub fn transmittance(
        self,
        wavelength_nm: f64,
        distance_m: f64,
    ) -> Result<f64, DielectricError> {
        if !distance_m.is_finite() || distance_m < 0.0 {
            return Err(DielectricError::InvalidAbsorption);
        }
        let extinction = self.extinction_per_m(wavelength_nm)?;
        Ok(det::exp(-extinction * distance_m))
    }

    /// Canonical source parameters for content identity and reporting.
    #[must_use]
    pub const fn parameters(self) -> BeerLambertParameters {
        match self.model {
            AbsorptionModel::Clear => BeerLambertParameters::Clear,
            AbsorptionModel::Constant(extinction_per_m) => {
                BeerLambertParameters::Constant { extinction_per_m }
            }
            AbsorptionModel::ReferenceRgb {
                linear_rgb,
                distance_m,
                ..
            } => BeerLambertParameters::ReferenceRgb {
                linear_rgb,
                distance_m,
            },
        }
    }
}

/// Provenance class attached to a glass definition.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum GlassProvenance {
    /// Caller-supplied parameters; their external source is not asserted here.
    Custom,
    /// Representative borosilicate-like visual preset v1, not measured stock.
    RepresentativeBorosilicateV1,
    /// Representative crown/low-iron-like visual preset v1, not measured stock.
    RepresentativeCrownV1,
}

/// A validated homogeneous dielectric medium.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DielectricGlass {
    ior: CauchyIor,
    absorption: BeerLambertAbsorption,
    provenance: GlassProvenance,
}

impl DielectricGlass {
    /// Construct a custom homogeneous medium from already validated parts.
    #[must_use]
    pub const fn new(
        ior: CauchyIor,
        absorption: BeerLambertAbsorption,
        provenance: GlassProvenance,
    ) -> Self {
        Self {
            ior,
            absorption,
            provenance,
        }
    }

    /// Representative borosilicate-like product glass.  The Cauchy fit and
    /// mild tint are visual starting values, not a measured material claim.
    #[must_use]
    pub fn representative_borosilicate() -> Self {
        let linear_rgb = [0.990, 0.996, 0.992];
        Self::new(
            CauchyIor {
                a: 1.4580,
                b_um2: 0.003_54,
                c_um4: 0.0,
            },
            BeerLambertAbsorption {
                model: AbsorptionModel::ReferenceRgb {
                    spectrum: lift_rgb(linear_rgb),
                    linear_rgb,
                    distance_m: 0.010,
                },
            },
            GlassProvenance::RepresentativeBorosilicateV1,
        )
    }

    /// Representative crown/low-iron-like product glass.  This is a visual
    /// preset with explicit provenance, not catalog or calibration data.
    #[must_use]
    pub fn representative_crown() -> Self {
        let linear_rgb = [0.986, 0.997, 0.992];
        Self::new(
            CauchyIor {
                a: 1.5046,
                b_um2: 0.004_20,
                c_um4: 0.0,
            },
            BeerLambertAbsorption {
                model: AbsorptionModel::ReferenceRgb {
                    spectrum: lift_rgb(linear_rgb),
                    linear_rgb,
                    distance_m: 0.010,
                },
            },
            GlassProvenance::RepresentativeCrownV1,
        )
    }

    /// Spectral phase-index law.
    #[must_use]
    pub const fn ior(self) -> CauchyIor {
        self.ior
    }

    /// Homogeneous absorption law.
    #[must_use]
    pub const fn absorption(self) -> BeerLambertAbsorption {
        self.absorption
    }

    /// Declared provenance class.
    #[must_use]
    pub const fn provenance(self) -> GlassProvenance {
        self.provenance
    }
}

/// Smooth or isotropic-GGX dielectric boundary treatment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DielectricSurface {
    roughness_alpha: Option<f64>,
}

impl DielectricSurface {
    /// Ideal smooth (delta) dielectric interface.
    pub const SMOOTH: Self = Self {
        roughness_alpha: None,
    };

    /// Polished-product-glass starting point (`alpha = 0.06`). This is a
    /// look-development convenience, not a measured surface claim.
    pub const POLISHED: Self = Self {
        roughness_alpha: Some(0.06),
    };

    /// Construct an isotropic GGX dielectric interface.
    pub fn try_rough(roughness_alpha: f64) -> Result<Self, DielectricError> {
        validate_roughness(roughness_alpha)?;
        Ok(Self {
            roughness_alpha: Some(roughness_alpha),
        })
    }

    /// `None` for a smooth delta interface, otherwise the GGX alpha.
    #[must_use]
    pub const fn roughness_alpha(self) -> Option<f64> {
        self.roughness_alpha
    }

    /// Whether the interface is a smooth delta distribution.
    #[must_use]
    pub const fn is_delta(self) -> bool {
        self.roughness_alpha.is_none()
    }
}

/// Exact unpolarized Fresnel result for one dielectric interface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FresnelDielectric {
    /// Power reflectance in `[0, 1]`.
    pub reflectance: f64,
    /// Positive transmitted-angle cosine, or zero under total internal
    /// reflection.
    pub transmitted_cosine: f64,
    /// Whether Snell's law admits no propagating transmitted ray.
    pub total_internal_reflection: bool,
}

/// Evaluate exact unpolarized dielectric Fresnel for an incident cosine in
/// `[0, 1]` and admitted absolute indices `eta_i`, `eta_t` in `[1, 3]`.
pub fn fresnel_dielectric(
    incident_cosine: f64,
    eta_i: f64,
    eta_t: f64,
) -> Result<FresnelDielectric, DielectricError> {
    validate_interface_scalar(incident_cosine, eta_i, eta_t)?;
    if incident_cosine < 0.0 || incident_cosine > 1.0 {
        return Err(DielectricError::InvalidInterface);
    }
    if indices_exactly_equal(eta_i, eta_t) {
        return Ok(FresnelDielectric {
            reflectance: 0.0,
            transmitted_cosine: incident_cosine,
            total_internal_reflection: false,
        });
    }
    let eta = eta_i / eta_t;
    let sin2_t = eta * eta * (1.0 - incident_cosine * incident_cosine).max(0.0);
    if sin2_t >= 1.0 {
        return Ok(FresnelDielectric {
            reflectance: 1.0,
            transmitted_cosine: 0.0,
            total_internal_reflection: true,
        });
    }
    let transmitted_cosine = (1.0 - sin2_t).max(0.0).sqrt();
    let parallel_numerator = eta_t * incident_cosine - eta_i * transmitted_cosine;
    let parallel_denominator = eta_t * incident_cosine + eta_i * transmitted_cosine;
    let perpendicular_numerator = eta_i * incident_cosine - eta_t * transmitted_cosine;
    let perpendicular_denominator = eta_i * incident_cosine + eta_t * transmitted_cosine;
    if parallel_denominator <= 0.0 || perpendicular_denominator <= 0.0 {
        return Err(DielectricError::InvalidInterface);
    }
    let parallel = parallel_numerator / parallel_denominator;
    let perpendicular = perpendicular_numerator / perpendicular_denominator;
    // Both ratios are bounded in [-1, 1], so the textbook half-sum cannot
    // overflow. Preserve this explicit Fresnel form rather than changing the
    // new transport bit semantics to satisfy the generic midpoint lint.
    #[allow(clippy::manual_midpoint)]
    let reflectance = 0.5 * (parallel * parallel + perpendicular * perpendicular);
    Ok(FresnelDielectric {
        reflectance: reflectance.clamp(0.0, 1.0),
        transmitted_cosine,
        total_internal_reflection: false,
    })
}

/// The sampled side of a dielectric interface.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DielectricEvent {
    /// Reflection remains in the incident medium.
    Reflection,
    /// Transmission crosses into the target medium.
    Transmission,
}

/// One smooth dielectric sample at one wavelength.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmoothDielectricSample {
    /// Unit outgoing direction, pointing away from the interface.
    pub direction: Vec3,
    /// Chosen reflection or transmission event.
    pub event: DielectricEvent,
    /// Discrete event probability under the supplied interface.
    pub probability: f64,
    /// Radiance-transport multiplier after division by event probability.
    pub radiance_weight: f64,
}

/// Sample an ideal smooth dielectric. `normal` and `wo` are unit vectors,
/// `normal.dot(wo) > 0`, and both directions point away from the surface.
pub fn sample_smooth_dielectric(
    normal: Vec3,
    wo: Vec3,
    eta_i: f64,
    eta_t: f64,
    event_sample: f64,
) -> Result<SmoothDielectricSample, DielectricError> {
    validate_direction_frame(normal, wo)?;
    validate_interface_scalar(normal.dot(wo), eta_i, eta_t)?;
    if !event_sample.is_finite() || !(0.0..1.0).contains(&event_sample) {
        return Err(DielectricError::InvalidInterface);
    }
    let normal = normalize_admitted(normal);
    let wo = normalize_admitted(wo);
    if indices_exactly_equal(eta_i, eta_t) {
        return Ok(SmoothDielectricSample {
            direction: wo.scale(-1.0),
            event: DielectricEvent::Transmission,
            probability: 1.0,
            radiance_weight: 1.0,
        });
    }
    // Normalizing two binary64 vectors does not guarantee their recomputed dot
    // is <= 1 by the last ulp. The frame validator already established the
    // geometric domain, so canonicalize that roundoff before Fresnel.
    let incident_cosine = normal.dot(wo).clamp(0.0, 1.0);
    let fresnel = fresnel_dielectric(incident_cosine, eta_i, eta_t)?;
    if fresnel.total_internal_reflection || event_sample < fresnel.reflectance {
        let direction = reflect(wo, normal);
        return Ok(SmoothDielectricSample {
            direction,
            event: DielectricEvent::Reflection,
            probability: fresnel.reflectance,
            radiance_weight: 1.0,
        });
    }
    let direction = refract(wo, normal, eta_i, eta_t, fresnel.transmitted_cosine)?;
    let eta_ratio = eta_i / eta_t;
    Ok(SmoothDielectricSample {
        direction,
        event: DielectricEvent::Transmission,
        probability: 1.0 - fresnel.reflectance,
        radiance_weight: eta_ratio * eta_ratio,
    })
}

/// Continuous rough-dielectric value and the matching visible-normal sampling
/// density at one wavelength. Smooth equal-IOR null boundaries are delta
/// distributions and therefore evaluate to zero in solid angle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RoughDielectricEvaluation {
    /// Radiance-mode BSDF value.
    pub value: f64,
    /// Solid-angle density of the view-conditioned GGX visible-normal plus
    /// Fresnel branch sampler.
    pub pdf: f64,
    /// Reflection or transmission hemisphere.
    pub event: DielectricEvent,
}

/// One rough dielectric sample at one wavelength.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RoughDielectricSample {
    /// Unit outgoing direction, pointing away from the interface.
    pub direction: Vec3,
    /// Reflection or transmission event.
    pub event: DielectricEvent,
    /// Matching solid-angle density for a rough event. Equal-IOR null
    /// transmission stores its discrete probability `1` here and sets
    /// [`Self::delta`], so callers must never feed that value to solid-angle
    /// MIS.
    pub pdf: f64,
    /// BSDF value at the sampled direction.
    pub value: f64,
    /// `value * abs(cos_theta_i) / pdf` in radiance transport.
    pub radiance_weight: f64,
    /// True only for the equal-IOR null-interface special case.
    pub delta: bool,
}

/// Evaluate an isotropic GGX dielectric BSDF and its matching sampler density.
pub fn evaluate_rough_dielectric(
    normal: Vec3,
    wo: Vec3,
    wi: Vec3,
    eta_i: f64,
    eta_t: f64,
    roughness_alpha: f64,
) -> Result<RoughDielectricEvaluation, DielectricError> {
    validate_direction_frame(normal, wo)?;
    validate_unit(wi)?;
    validate_interface_scalar(normal.dot(wo), eta_i, eta_t)?;
    validate_roughness(roughness_alpha)?;
    let normal = normalize_admitted(normal);
    let wo = normalize_admitted(wo);
    let wi = normalize_admitted(wi);
    let cos_o = normal.dot(wo).clamp(0.0, 1.0);
    let cos_i = normal.dot(wi).clamp(-1.0, 1.0);
    if cos_i == 0.0 {
        return Ok(RoughDielectricEvaluation {
            value: 0.0,
            pdf: 0.0,
            event: DielectricEvent::Transmission,
        });
    }
    let reflection = cos_i > 0.0;
    let event = if reflection {
        DielectricEvent::Reflection
    } else {
        DielectricEvent::Transmission
    };
    if indices_exactly_equal(eta_i, eta_t) && !reflection {
        return Ok(RoughDielectricEvaluation {
            value: 0.0,
            pdf: 0.0,
            event,
        });
    }

    let eta = eta_t / eta_i;
    let half_sum = if reflection {
        add(wo, wi)
    } else {
        add(wo, wi.scale(eta))
    };
    let half_norm = scale_safe_norm(half_sum);
    if !half_norm.is_finite() || half_norm <= 0.0 {
        return Ok(RoughDielectricEvaluation {
            value: 0.0,
            pdf: 0.0,
            event,
        });
    }
    let mut micro_normal = half_sum.scale(1.0 / half_norm);
    if normal.dot(micro_normal) < 0.0 {
        micro_normal = micro_normal.scale(-1.0);
    }
    let n_dot_m = normal.dot(micro_normal).clamp(-1.0, 1.0);
    let wo_dot_m = wo.dot(micro_normal).clamp(-1.0, 1.0);
    let wi_dot_m = wi.dot(micro_normal).clamp(-1.0, 1.0);
    if n_dot_m <= 0.0 || wo_dot_m <= 0.0 {
        return Ok(RoughDielectricEvaluation {
            value: 0.0,
            pdf: 0.0,
            event,
        });
    }
    let fresnel = fresnel_dielectric(wo_dot_m.clamp(0.0, 1.0), eta_i, eta_t)?;
    let distribution = ggx_d(roughness_alpha, n_dot_m);
    let geometry = smith_g1(roughness_alpha, cos_o) * smith_g1(roughness_alpha, cos_i.abs());

    if reflection {
        if wi_dot_m <= 0.0 {
            return Ok(RoughDielectricEvaluation {
                value: 0.0,
                pdf: 0.0,
                event,
            });
        }
        let value = distribution * geometry * fresnel.reflectance / (4.0 * cos_o * cos_i);
        let visible_normal_pdf = distribution * smith_g1(roughness_alpha, cos_o) * wo_dot_m / cos_o;
        let pdf = visible_normal_pdf * fresnel.reflectance / (4.0 * wo_dot_m);
        return finite_evaluation(value, pdf, event);
    }

    if wi_dot_m >= 0.0 || fresnel.total_internal_reflection {
        return Ok(RoughDielectricEvaluation {
            value: 0.0,
            pdf: 0.0,
            event,
        });
    }
    // `half_sum = wo + eta*wi = denominator*m`. Retain its scale before
    // normalization instead of subtracting the two nearly equal dot terms a
    // second time. This keeps distinct adjacent-f64 indices finite at normal
    // incidence while preserving the Walter Jacobian exactly in real
    // arithmetic.
    let denominator2 = half_norm * half_norm;
    if !denominator2.is_finite() || denominator2 <= 0.0 {
        return Ok(RoughDielectricEvaluation {
            value: 0.0,
            pdf: 0.0,
            event,
        });
    }
    // Walter et al. microfacet transmission in radiance mode.  Omitting the
    // eta^2 numerator from `value` is the radiance Jacobian; the sampler PDF
    // retains eta^2, so a smooth-limit slab gains eta_i^2/eta_t^2 on entry and
    // the reciprocal on exit.
    let value = (1.0 - fresnel.reflectance) * distribution * geometry * (wi_dot_m * wo_dot_m).abs()
        / (cos_o * cos_i.abs() * denominator2);
    let jacobian = eta * eta * wi_dot_m.abs() / denominator2;
    let visible_normal_pdf = distribution * smith_g1(roughness_alpha, cos_o) * wo_dot_m / cos_o;
    let pdf = visible_normal_pdf * (1.0 - fresnel.reflectance) * jacobian;
    finite_evaluation(value, pdf, event)
}

/// Sample an isotropic GGX dielectric using Heitz's view-conditioned visible
/// normal distribution and an exact Fresnel branch decision. A valid but
/// zero-contribution microfacet draw returns `Ok(None)`.
#[allow(clippy::too_many_arguments)]
pub fn sample_rough_dielectric(
    normal: Vec3,
    wo: Vec3,
    eta_i: f64,
    eta_t: f64,
    roughness_alpha: f64,
    microfacet_u: f64,
    microfacet_v: f64,
    event_sample: f64,
) -> Result<Option<RoughDielectricSample>, DielectricError> {
    validate_direction_frame(normal, wo)?;
    validate_interface_scalar(normal.dot(wo), eta_i, eta_t)?;
    validate_roughness(roughness_alpha)?;
    for sample in [microfacet_u, microfacet_v, event_sample] {
        if !sample.is_finite() || !(0.0..1.0).contains(&sample) {
            return Err(DielectricError::InvalidInterface);
        }
    }
    let normal = normalize_admitted(normal);
    let wo = normalize_admitted(wo);
    if indices_exactly_equal(eta_i, eta_t) {
        return Ok(Some(RoughDielectricSample {
            direction: wo.scale(-1.0),
            event: DielectricEvent::Transmission,
            pdf: 1.0,
            value: 0.0,
            radiance_weight: 1.0,
            delta: true,
        }));
    }

    let micro_normal =
        sample_ggx_visible_normal(normal, wo, roughness_alpha, microfacet_u, microfacet_v)?;
    let wo_dot_m = wo.dot(micro_normal);
    if wo_dot_m <= 0.0 {
        return Ok(None);
    }
    let fresnel = fresnel_dielectric(wo_dot_m.clamp(0.0, 1.0), eta_i, eta_t)?;
    let (direction, event) =
        if fresnel.total_internal_reflection || event_sample < fresnel.reflectance {
            (reflect(wo, micro_normal), DielectricEvent::Reflection)
        } else {
            (
                refract(wo, micro_normal, eta_i, eta_t, fresnel.transmitted_cosine)?,
                DielectricEvent::Transmission,
            )
        };
    let cos_i = normal.dot(direction);
    if (event == DielectricEvent::Reflection && cos_i <= 0.0)
        || (event == DielectricEvent::Transmission && cos_i >= 0.0)
    {
        return Ok(None);
    }
    let evaluation =
        evaluate_rough_dielectric(normal, wo, direction, eta_i, eta_t, roughness_alpha)?;
    if evaluation.pdf <= 0.0 || evaluation.value < 0.0 {
        return Ok(None);
    }
    let radiance_weight = evaluation.value * cos_i.abs() / evaluation.pdf;
    if !radiance_weight.is_finite() || radiance_weight < 0.0 {
        return Err(DielectricError::InvalidInterface);
    }
    Ok(Some(RoughDielectricSample {
        direction,
        event,
        pdf: evaluation.pdf,
        value: evaluation.value,
        radiance_weight,
        delta: false,
    }))
}

fn validate_wavelength(wavelength_nm: f64) -> Result<(), DielectricError> {
    if wavelength_nm.is_finite() && (LAMBDA_MIN..=LAMBDA_MAX).contains(&wavelength_nm) {
        Ok(())
    } else {
        Err(DielectricError::InvalidWavelength)
    }
}

fn validate_interface_scalar(
    incident_cosine: f64,
    eta_i: f64,
    eta_t: f64,
) -> Result<(), DielectricError> {
    if incident_cosine.is_finite()
        && eta_i.is_finite()
        && eta_t.is_finite()
        && (1.0..=MAX_GLASS_IOR).contains(&eta_i)
        && (1.0..=MAX_GLASS_IOR).contains(&eta_t)
    {
        Ok(())
    } else {
        Err(DielectricError::InvalidInterface)
    }
}

fn validate_roughness(alpha: f64) -> Result<(), DielectricError> {
    if alpha.is_finite() && (MIN_ROUGHNESS_ALPHA..=MAX_ROUGHNESS_ALPHA).contains(&alpha) {
        Ok(())
    } else {
        Err(DielectricError::InvalidRoughness)
    }
}

fn validate_direction_frame(normal: Vec3, wo: Vec3) -> Result<(), DielectricError> {
    validate_unit(normal)?;
    validate_unit(wo)?;
    let cosine = normal.dot(wo);
    if cosine.is_finite() && cosine > 0.0 && cosine <= 1.0 + UNIT_TOLERANCE {
        Ok(())
    } else {
        Err(DielectricError::InvalidDirection)
    }
}

fn validate_unit(direction: Vec3) -> Result<(), DielectricError> {
    if !direction.x.is_finite() || !direction.y.is_finite() || !direction.z.is_finite() {
        return Err(DielectricError::InvalidDirection);
    }
    let norm = direction.norm();
    if norm.is_finite() && (norm - 1.0).abs() <= UNIT_TOLERANCE {
        Ok(())
    } else {
        Err(DielectricError::InvalidDirection)
    }
}

fn normalize_admitted(direction: Vec3) -> Vec3 {
    direction.scale(1.0 / direction.norm())
}

fn indices_exactly_equal(eta_i: f64, eta_t: f64) -> bool {
    eta_i.to_bits() == eta_t.to_bits()
}

fn scale_safe_norm(vector: Vec3) -> f64 {
    let scale = vector.x.abs().max(vector.y.abs()).max(vector.z.abs());
    if scale == 0.0 {
        0.0
    } else if !scale.is_finite() {
        f64::INFINITY
    } else {
        let scaled = vector.scale(1.0 / scale);
        scale * scaled.dot(scaled).sqrt()
    }
}

fn reflect(wo: Vec3, normal: Vec3) -> Vec3 {
    let cosine = wo.dot(normal);
    add(normal.scale(2.0 * cosine), wo.scale(-1.0))
}

fn refract(
    wo: Vec3,
    normal: Vec3,
    eta_i: f64,
    eta_t: f64,
    transmitted_cosine: f64,
) -> Result<Vec3, DielectricError> {
    let incident_cosine = wo.dot(normal);
    let eta = eta_i / eta_t;
    let direction = add(
        normal.scale(eta * incident_cosine - transmitted_cosine),
        wo.scale(-eta),
    );
    let norm = direction.norm();
    if !norm.is_finite() || norm <= 0.0 {
        return Err(DielectricError::InvalidDirection);
    }
    let unit = direction.scale(1.0 / norm);
    if normal.dot(unit) >= 0.0 {
        Err(DielectricError::InvalidDirection)
    } else {
        Ok(unit)
    }
}

fn add(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}

fn finite_evaluation(
    value: f64,
    pdf: f64,
    event: DielectricEvent,
) -> Result<RoughDielectricEvaluation, DielectricError> {
    if value.is_finite() && pdf.is_finite() && value >= 0.0 && pdf >= 0.0 {
        Ok(RoughDielectricEvaluation { value, pdf, event })
    } else {
        Err(DielectricError::InvalidInterface)
    }
}

fn ggx_d(alpha: f64, cos_m: f64) -> f64 {
    if cos_m <= 0.0 {
        return 0.0;
    }
    let a2 = alpha * alpha;
    let c2 = cos_m * cos_m;
    let denominator = c2 * (a2 - 1.0) + 1.0;
    a2 / (PI * denominator * denominator)
}

fn smith_g1(alpha: f64, cosine: f64) -> f64 {
    if cosine <= 0.0 {
        return 0.0;
    }
    let a2 = alpha * alpha;
    2.0 * cosine / (cosine + (a2 + (1.0 - a2) * cosine * cosine).sqrt())
}

fn basis(normal: Vec3) -> (Vec3, Vec3) {
    // Duff et al.'s cancellation-free all-sphere ONB. The older Frisvad
    // south-pole shortcut was exact only at (0, 0, -1); applying it to a cap
    // produced non-unit microfacet normals and a sampling/PDF mismatch.
    let sign = if normal.z < 0.0 { -1.0 } else { 1.0 };
    let a = -1.0 / (sign + normal.z);
    let b = normal.x * normal.y * a;
    (
        Vec3::new(
            1.0 + sign * normal.x * normal.x * a,
            sign * b,
            -sign * normal.x,
        ),
        Vec3::new(b, sign + normal.y * normal.y * a, -normal.y),
    )
}

fn to_world(normal: Vec3, local: [f64; 3]) -> Vec3 {
    let (tangent, bitangent) = basis(normal);
    Vec3::new(
        tangent.x * local[0] + bitangent.x * local[1] + normal.x * local[2],
        tangent.y * local[0] + bitangent.y * local[1] + normal.y * local[2],
        tangent.z * local[0] + bitangent.z * local[1] + normal.z * local[2],
    )
}

/// Heitz's isotropic GGX visible-normal warp (JCGT 2018). Conditioning the
/// microfacet proposal on `wo` avoids spending most grazing-incidence draws on
/// facets hidden from the incident direction. The returned density is
/// evaluated analytically in [`evaluate_rough_dielectric`].
fn sample_ggx_visible_normal(
    normal: Vec3,
    wo: Vec3,
    alpha: f64,
    u1: f64,
    u2: f64,
) -> Result<Vec3, DielectricError> {
    let (tangent, bitangent) = basis(normal);
    let local_wo = [tangent.dot(wo), bitangent.dot(wo), normal.dot(wo)];
    let stretched = [alpha * local_wo[0], alpha * local_wo[1], local_wo[2]];
    let stretched_norm =
        (stretched[0] * stretched[0] + stretched[1] * stretched[1] + stretched[2] * stretched[2])
            .sqrt();
    if !stretched_norm.is_finite() || stretched_norm <= 0.0 {
        return Err(DielectricError::InvalidInterface);
    }
    let view = stretched.map(|component| component / stretched_norm);

    let lensq = view[0] * view[0] + view[1] * view[1];
    let frame_1 = if lensq > 0.0 {
        let inverse_len = 1.0 / lensq.sqrt();
        [-view[1] * inverse_len, view[0] * inverse_len, 0.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    let frame_2 = [
        view[1] * frame_1[2] - view[2] * frame_1[1],
        view[2] * frame_1[0] - view[0] * frame_1[2],
        view[0] * frame_1[1] - view[1] * frame_1[0],
    ];

    let radius = u1.sqrt();
    let phi = 2.0 * PI * u2;
    let disk_1 = radius * det::cos(phi);
    let mut disk_2 = radius * det::sin(phi);
    let blend = 0.5 * (1.0 + view[2]);
    disk_2 = (1.0 - blend) * (1.0 - disk_1 * disk_1).max(0.0).sqrt() + blend * disk_2;
    let projected = (1.0 - disk_1 * disk_1 - disk_2 * disk_2).max(0.0).sqrt();
    let stretched_normal = [
        disk_1 * frame_1[0] + disk_2 * frame_2[0] + projected * view[0],
        disk_1 * frame_1[1] + disk_2 * frame_2[1] + projected * view[1],
        disk_1 * frame_1[2] + disk_2 * frame_2[2] + projected * view[2],
    ];

    let unstretched = [
        alpha * stretched_normal[0],
        alpha * stretched_normal[1],
        stretched_normal[2].max(0.0),
    ];
    let normal_length = (unstretched[0] * unstretched[0]
        + unstretched[1] * unstretched[1]
        + unstretched[2] * unstretched[2])
        .sqrt();
    if !normal_length.is_finite() || normal_length <= 0.0 {
        return Err(DielectricError::InvalidInterface);
    }
    let local_normal = unstretched.map(|component| component / normal_length);
    let micro_normal = to_world(normal, local_normal);
    if !micro_normal.x.is_finite()
        || !micro_normal.y.is_finite()
        || !micro_normal.z.is_finite()
        || wo.dot(micro_normal) <= 0.0
    {
        return Err(DielectricError::InvalidInterface);
    }
    Ok(micro_normal)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(observed: f64, expected: f64, tolerance: f64, context: &str) {
        assert!(
            (observed - expected).abs() <= tolerance,
            "{context}: observed {observed:.17e}, expected {expected:.17e}"
        );
    }

    #[test]
    fn normal_incidence_fresnel_matches_closed_form() {
        let result = fresnel_dielectric(1.0, 1.0, 1.5).expect("valid interface");
        assert_close(result.reflectance, 0.04, 2.0e-15, "normal Fresnel");
        assert_close(result.transmitted_cosine, 1.0, 0.0, "normal Snell");
        assert!(!result.total_internal_reflection);
    }

    #[test]
    fn fresnel_angle_table_matches_independent_known_answers() {
        for (degrees, expected_cosine, expected_reflectance) in [
            (30.0_f64, 0.942_809_041_582_063_4, 0.041_522_625_975_821_52),
            (60.0_f64, 0.816_496_580_927_726, 0.089_186_712_802_212_74),
            (89.9_f64, 0.745_356_900_692_131_5, 0.989_911_876_919_429_7),
        ] {
            let incident_cosine = degrees.to_radians().cos();
            let result = fresnel_dielectric(incident_cosine, 1.0, 1.5).expect("air/glass");
            assert_close(
                result.transmitted_cosine,
                expected_cosine,
                3.0e-15,
                "tabulated Snell cosine",
            );
            assert_close(
                result.reflectance,
                expected_reflectance,
                3.0e-15,
                "tabulated Fresnel",
            );
        }
        let equal_grazing = fresnel_dielectric(0.0, 1.5, 1.5).expect("equal-IOR grazing");
        assert!(!equal_grazing.total_internal_reflection);
        assert_eq!(equal_grazing.reflectance.to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn snell_and_total_internal_reflection_are_exactly_classified() {
        let theta_i = core::f64::consts::FRAC_PI_6;
        let incident_cosine = theta_i.cos();
        let result = fresnel_dielectric(incident_cosine, 1.0, 1.5).expect("air/glass");
        let expected_cosine = (1.0 - (theta_i.sin() / 1.5).powi(2)).sqrt();
        assert_close(
            result.transmitted_cosine,
            expected_cosine,
            3.0e-15,
            "Snell cosine",
        );

        let above_critical = fresnel_dielectric(0.5, 1.5, 1.0).expect("glass/air");
        assert!(above_critical.total_internal_reflection);
        assert_eq!(above_critical.reflectance.to_bits(), 1.0_f64.to_bits());

        let critical_cosine = (1.0_f64 - (1.0_f64 / 1.5).powi(2)).sqrt();
        let at_critical =
            fresnel_dielectric(critical_cosine, 1.5, 1.0).expect("critical glass/air");
        assert!(at_critical.total_internal_reflection);
        assert_eq!(at_critical.reflectance.to_bits(), 1.0_f64.to_bits());
        let just_above_angle =
            fresnel_dielectric(f64::from_bits(critical_cosine.to_bits() - 1), 1.5, 1.0)
                .expect("one ulp above critical angle");
        assert!(just_above_angle.total_internal_reflection);
        let just_below_angle =
            fresnel_dielectric(f64::from_bits(critical_cosine.to_bits() + 1), 1.5, 1.0)
                .expect("one ulp below critical angle");
        assert!(!just_below_angle.total_internal_reflection);
        assert!(just_below_angle.transmitted_cosine > 0.0);
    }

    #[test]
    fn smooth_refraction_obeys_vector_snell_and_equal_ior_is_null() {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let theta_i = core::f64::consts::FRAC_PI_6;
        let wo = Vec3::new(theta_i.sin(), 0.0, theta_i.cos());
        let sample =
            sample_smooth_dielectric(normal, wo, 1.0, 1.5, 0.99).expect("transmission sample");
        assert_eq!(sample.event, DielectricEvent::Transmission);
        assert_close(
            sample.direction.x,
            -theta_i.sin() / 1.5,
            3.0e-15,
            "vector Snell",
        );
        assert!(sample.direction.z < 0.0);
        assert_close(sample.radiance_weight, 1.0 / 2.25, 2.0e-15, "eta scale");

        let null =
            sample_smooth_dielectric(normal, wo, 1.0, 1.0, 0.0).expect("equal-IOR null boundary");
        assert_eq!(null.event, DielectricEvent::Transmission);
        assert_eq!(null.direction, normalize_admitted(wo).scale(-1.0));
        assert_eq!(null.probability.to_bits(), 1.0_f64.to_bits());
        assert_eq!(null.radiance_weight.to_bits(), 1.0_f64.to_bits());
    }

    #[test]
    fn lossless_smooth_slab_cancels_radiance_eta_factors() {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let wo = normal;
        let entry =
            sample_smooth_dielectric(normal, wo, 1.0, 1.5, 0.5).expect("air-to-glass transmission");
        let exit =
            sample_smooth_dielectric(normal, wo, 1.5, 1.0, 0.5).expect("glass-to-air transmission");
        assert_eq!(entry.event, DielectricEvent::Transmission);
        assert_eq!(exit.event, DielectricEvent::Transmission);
        assert_close(
            entry.radiance_weight * exit.radiance_weight,
            1.0,
            2.0e-15,
            "lossless slab eta cancellation",
        );
    }

    #[test]
    fn beer_lambert_path_scaling_is_multiplicative() {
        let absorption = BeerLambertAbsorption::try_constant(17.0).expect("constant extinction");
        let one = absorption.transmittance(550.0, 0.012).expect("one length");
        let two = absorption.transmittance(550.0, 0.024).expect("two lengths");
        assert_close(two, one * one, 4.0e-15, "Beer-Lambert doubling");
        assert_eq!(
            absorption
                .transmittance(550.0, 0.0)
                .expect("zero length")
                .to_bits(),
            1.0_f64.to_bits()
        );
    }

    #[test]
    fn cauchy_and_absorption_admission_refuse_nonphysical_inputs() {
        assert_eq!(
            CauchyIor::try_new(0.99, 0.0, 0.0),
            Err(DielectricError::InvalidIor)
        );
        assert_eq!(
            CauchyIor::try_new(1.5, -0.01, 0.0),
            Err(DielectricError::InvalidIor)
        );
        assert_eq!(
            BeerLambertAbsorption::try_constant(-1.0),
            Err(DielectricError::InvalidAbsorption)
        );
        assert_eq!(
            BeerLambertAbsorption::try_from_rgb_transmittance([1.1, 0.9, 0.9], 0.01),
            Err(DielectricError::InvalidAbsorption)
        );
        assert_eq!(
            DielectricSurface::try_rough(0.0),
            Err(DielectricError::InvalidRoughness)
        );
        for invalid_eta in [0.0, -1.0, 0.99, 3.01, f64::MAX, f64::INFINITY] {
            assert_eq!(
                sample_smooth_dielectric(
                    Vec3::new(0.0, 0.0, 1.0),
                    Vec3::new(0.0, 0.0, 1.0),
                    invalid_eta,
                    invalid_eta,
                    0.5,
                ),
                Err(DielectricError::InvalidInterface)
            );
        }
        for boundary_eta in [1.0, MAX_GLASS_IOR] {
            assert!(fresnel_dielectric(1.0, boundary_eta, boundary_eta).is_ok());
        }
        let signed_zero = CauchyIor::try_new(1.5, -0.0, -0.0).expect("signed zero");
        assert!(!signed_zero.is_dispersive());
        assert_eq!(signed_zero.coefficients()[1].to_bits(), 0.0_f64.to_bits());
        assert_eq!(
            BeerLambertAbsorption::try_constant(-0.0)
                .expect("signed zero extinction")
                .parameters(),
            BeerLambertParameters::Clear
        );
    }

    #[test]
    fn admitted_non_axis_directions_are_normalized_before_every_sampling_return() {
        // A second normalization of this ordinary non-axis vector self-dots to
        // 1 + 1 ulp on binary64, so it exercises both direction normalization
        // and the post-validation Fresnel cosine clamp.
        let near_unit = Vec3::new(
            0.792_305_471_084_122_2,
            0.238_331_816_145_395_78,
            0.561_649_344_255_830_4,
        );
        let smooth = sample_smooth_dielectric(near_unit, near_unit, 1.0, 1.5, 0.5)
            .expect("admitted near-unit smooth frame");
        assert_close(
            smooth.direction.norm(),
            1.0,
            2.0e-15,
            "normalized smooth direction",
        );
        let rough = sample_rough_dielectric(near_unit, near_unit, 1.0, 1.5, 0.2, 0.37, 0.61, 0.0)
            .expect("admitted near-unit rough frame")
            .expect("rough reflection");
        assert_close(
            rough.direction.norm(),
            1.0,
            2.0e-15,
            "normalized rough direction",
        );

        let smooth_null = sample_smooth_dielectric(near_unit, near_unit, 1.5, 1.5, 0.0)
            .expect("equal-IOR smooth early return");
        assert_eq!(smooth_null.event, DielectricEvent::Transmission);
        assert_close(
            smooth_null.direction.norm(),
            1.0,
            2.0e-15,
            "normalized smooth null direction",
        );
        let rough_null =
            sample_rough_dielectric(near_unit, near_unit, 1.5, 1.5, 0.2, 0.37, 0.61, 0.0)
                .expect("equal-IOR rough early return")
                .expect("equal-IOR rough null sample");
        assert!(rough_null.delta);
        assert_close(
            rough_null.direction.norm(),
            1.0,
            2.0e-15,
            "normalized rough null direction",
        );
    }

    #[test]
    fn adjacent_distinct_indices_retain_fresnel_and_stable_rough_transmission() {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let eta = 1.5_f64;
        let lower = f64::from_bits(eta.to_bits() - 1);
        let upper = f64::from_bits(eta.to_bits() + 1);
        assert!(indices_exactly_equal(eta, eta));
        for adjacent_eta in [lower, upper] {
            assert!(!indices_exactly_equal(eta, adjacent_eta));
            let grazing =
                fresnel_dielectric(0.0, eta, adjacent_eta).expect("adjacent-IOR grazing Fresnel");
            assert_eq!(grazing.reflectance.to_bits(), 1.0_f64.to_bits());
            assert_eq!(
                grazing.total_internal_reflection,
                eta > adjacent_eta,
                "adjacent-IOR grazing TIR classification"
            );
            let normal_fresnel =
                fresnel_dielectric(1.0, eta, adjacent_eta).expect("adjacent-IOR normal Fresnel");
            let expected_normal_reflectance = ((eta - adjacent_eta) / (eta + adjacent_eta)).powi(2);
            assert_close(
                normal_fresnel.reflectance,
                expected_normal_reflectance,
                8.0 * f64::EPSILON * expected_normal_reflectance,
                "adjacent-IOR normal Fresnel closed form",
            );

            let sample =
                sample_rough_dielectric(normal, normal, eta, adjacent_eta, 0.2, 0.0, 0.0, 0.5)
                    .expect("adjacent-IOR rough evaluation")
                    .expect("adjacent-IOR transmission remains contributing");
            assert!(!sample.delta);
            assert_eq!(sample.event, DielectricEvent::Transmission);
            assert!(sample.pdf.is_finite() && sample.pdf > 0.0);
            assert!(sample.value.is_finite() && sample.value > 0.0);
            assert!(sample.radiance_weight.is_finite() && sample.radiance_weight > 0.0);
            assert_close(
                sample.radiance_weight,
                (eta / adjacent_eta).powi(2),
                2.0e-15,
                "adjacent-IOR radiance factor",
            );
        }
    }

    #[test]
    fn representative_presets_are_dispersive_and_honestly_labeled() {
        let glass = DielectricGlass::representative_borosilicate();
        let blue = glass.ior().eval(420.0).expect("blue");
        let red = glass.ior().eval(700.0).expect("red");
        assert!(blue > red && red > 1.0);
        assert_eq!(
            glass.provenance(),
            GlassProvenance::RepresentativeBorosilicateV1
        );
        assert!(
            glass
                .absorption()
                .transmittance(550.0, 0.01)
                .expect("green transmission")
                > 0.9
        );
    }

    #[test]
    fn rough_samples_match_their_evaluation_and_radiance_jacobian() {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let wo = Vec3::new(0.3, 0.0, (1.0_f64 - 0.09).sqrt());
        for event_sample in [0.0, 0.9] {
            let sample =
                sample_rough_dielectric(normal, wo, 1.0, 1.5, 0.2, 0.17, 0.63, event_sample)
                    .expect("valid rough sample")
                    .expect("contributing microfacet");
            let evaluated = evaluate_rough_dielectric(normal, wo, sample.direction, 1.0, 1.5, 0.2)
                .expect("evaluate sample");
            assert_eq!(sample.event, evaluated.event);
            assert_close(sample.value, evaluated.value, 0.0, "rough value replay");
            assert_close(sample.pdf, evaluated.pdf, 0.0, "rough PDF replay");
            assert_close(
                sample.radiance_weight,
                evaluated.value * normal.dot(sample.direction).abs() / evaluated.pdf,
                2.0e-15,
                "rough throughput",
            );
        }
    }

    #[test]
    fn rough_visible_normal_pdf_mass_matches_sampler_acceptance() {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let cos_o = 0.2_f64;
        let wo = Vec3::new((1.0 - cos_o * cos_o).sqrt(), 0.0, cos_o);
        let alpha = 0.15;
        let sphere_rows = 128_u32;
        let sphere_columns = 256_u32;
        let mut pdf_sum = 0.0;
        for row in 0..sphere_rows {
            let z = -1.0 + 2.0 * (f64::from(row) + 0.5) / f64::from(sphere_rows);
            let radial = (1.0 - z * z).sqrt();
            for column in 0..sphere_columns {
                let phi = 2.0 * PI * (f64::from(column) + 0.5) / f64::from(sphere_columns);
                let wi = Vec3::new(radial * det::cos(phi), radial * det::sin(phi), z);
                pdf_sum += evaluate_rough_dielectric(normal, wo, wi, 1.0, 1.5, alpha)
                    .expect("finite directional evaluation")
                    .pdf;
            }
        }
        let sphere_count = f64::from(sphere_rows) * f64::from(sphere_columns);
        let integrated_mass = pdf_sum * 4.0 * PI / sphere_count;

        let side = 32_u32;
        let mut accepted = 0_u32;
        for microfacet_u in 0..side {
            for microfacet_v in 0..side {
                for event in 0..side {
                    accepted += u32::from(
                        sample_rough_dielectric(
                            normal,
                            wo,
                            1.0,
                            1.5,
                            alpha,
                            (f64::from(microfacet_u) + 0.5) / f64::from(side),
                            (f64::from(microfacet_v) + 0.5) / f64::from(side),
                            (f64::from(event) + 0.5) / f64::from(side),
                        )
                        .expect("finite visible-normal sample")
                        .is_some(),
                    );
                }
            }
        }
        let sampled_mass = f64::from(accepted) / f64::from(side * side * side);
        assert_close(
            integrated_mass,
            sampled_mass,
            1.5e-2,
            "rough VNDF directional PDF mass versus sampler acceptance",
        );
    }

    #[test]
    fn rough_normal_transmission_has_the_expected_eta_squared_limit_factor() {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let wo = normal;
        let sample = sample_rough_dielectric(normal, wo, 1.0, 1.5, 1.0e-4, 0.0, 0.0, 0.99)
            .expect("valid sample")
            .expect("transmission");
        assert_eq!(sample.event, DielectricEvent::Transmission);
        assert_close(
            sample.radiance_weight,
            1.0 / 2.25,
            3.0e-12,
            "rough smooth-limit eta factor",
        );
    }

    #[test]
    fn rough_walter_fixture_matches_independent_values() {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let wo = normal;
        let microfacet_u = 0.679_444_983_649_738_8;
        let transmission =
            sample_rough_dielectric(normal, wo, 1.0, 1.5, 0.25, microfacet_u, 0.0, 0.5)
                .expect("valid Walter fixture")
                .expect("transmission sample");
        assert_eq!(transmission.event, DielectricEvent::Transmission);
        assert_close(
            transmission.direction.x,
            -0.118_748_107_807_083_9,
            3.0e-15,
            "Walter transmission x",
        );
        assert_close(
            transmission.direction.z,
            -0.992_924_411_469_592_8,
            3.0e-15,
            "Walter transmission z",
        );
        assert_close(
            transmission.pdf,
            4.889_055_831_876_031,
            3.0e-14,
            "Walter transmission PDF",
        );
        assert_close(
            transmission.value,
            2.187_909_058_805_138,
            3.0e-14,
            "Walter transmission value",
        );
        assert_close(
            transmission.radiance_weight,
            0.444_345_163_824_704_57,
            3.0e-14,
            "Walter transmission weight",
        );

        let reflection =
            sample_rough_dielectric(normal, wo, 1.0, 1.5, 0.25, microfacet_u, 0.0, 0.0)
                .expect("valid Walter fixture")
                .expect("reflection sample");
        assert_eq!(reflection.event, DielectricEvent::Reflection);
        assert_close(
            reflection.direction.x,
            0.642_787_609_686_539_4,
            3.0e-15,
            "Walter reflection x",
        );
        assert_close(
            reflection.pdf,
            0.006_756_362_821_091_062_5,
            3.0e-16,
            "Walter reflection PDF",
        );
        assert_close(
            reflection.value,
            0.008_724_853_225_265_146,
            3.0e-16,
            "Walter reflection value",
        );
        assert_close(
            reflection.radiance_weight,
            0.989_234_223_683_053,
            3.0e-14,
            "Walter reflection weight",
        );
    }

    #[test]
    fn rough_sampling_frame_is_orthonormal_beside_both_poles() {
        let x = 4.0e-4_f64;
        let z = (1.0 - x * x).sqrt();
        let normals = [Vec3::new(x, 0.0, z), Vec3::new(x, 0.0, -z)];
        let mut samples = Vec::new();
        for normal in normals {
            let (tangent, bitangent) = basis(normal);
            assert_close(tangent.norm(), 1.0, 4.0e-15, "pole tangent norm");
            assert_close(bitangent.norm(), 1.0, 4.0e-15, "pole bitangent norm");
            assert_close(normal.dot(tangent), 0.0, 4.0e-15, "pole n dot t");
            assert_close(normal.dot(bitangent), 0.0, 4.0e-15, "pole n dot b");
            assert_close(tangent.dot(bitangent), 0.0, 4.0e-15, "pole t dot b");

            let sample = sample_rough_dielectric(normal, normal, 1.0, 1.5, 0.2, 0.37, 0.61, 0.0)
                .expect("pole-adjacent sample must retain a unit frame")
                .expect("pole-adjacent reflection contributes");
            assert_eq!(sample.event, DielectricEvent::Reflection);
            assert_close(
                sample.direction.norm(),
                1.0,
                4.0e-15,
                "sample direction norm",
            );
            let replay = evaluate_rough_dielectric(normal, normal, sample.direction, 1.0, 1.5, 0.2)
                .expect("pole-adjacent evaluation");
            assert_close(sample.value, replay.value, 0.0, "pole value replay");
            assert_close(sample.pdf, replay.pdf, 0.0, "pole PDF replay");
            samples.push(sample);
        }
        assert_close(
            samples[0].value,
            samples[1].value,
            2.0e-14,
            "north/south isotropic value",
        );
        assert_close(
            samples[0].pdf,
            samples[1].pdf,
            2.0e-14,
            "north/south isotropic PDF",
        );
    }

    #[test]
    fn rough_grazing_grid_stays_finite_and_nonnegative() {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        for cosine in [1.0, 0.5, 0.1, 0.01, 0.001] {
            let wo = Vec3::new((1.0_f64 - cosine * cosine).sqrt(), 0.0, cosine);
            for alpha in [1.0e-4, 0.01, 0.2, 1.0] {
                for event_sample in [0.0, 0.5, 0.999_999] {
                    if let Some(sample) = sample_rough_dielectric(
                        normal,
                        wo,
                        1.0,
                        1.5,
                        alpha,
                        0.37,
                        0.61,
                        event_sample,
                    )
                    .expect("finite grazing fixture")
                    {
                        assert!(sample.pdf.is_finite() && sample.pdf >= 0.0);
                        assert!(sample.value.is_finite() && sample.value >= 0.0);
                        assert!(
                            sample.radiance_weight.is_finite() && sample.radiance_weight >= 0.0
                        );
                    }
                }
            }
        }
    }
}

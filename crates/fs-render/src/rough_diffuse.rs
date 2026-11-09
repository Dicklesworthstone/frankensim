//! Reciprocal, energy-preserving Oren--Nayar (EON) diffuse reflection.
//!
//! This is the analytic EON model of Portsmouth, Kutz and Hill, JCGT 14(1),
//! equations 10--20: <https://jcgt.org/published/0014/01/06/>. The normalized
//! roughness is a model parameter in [0,1], not an RMS slope angle or a height
//! measurement. The spectral input is the single-scattering albedo of the
//! diffuse facets, not the resulting directional or hemispherical reflectance.
//! Multiple returns use the paper's reciprocal analytic approximation. This
//! does not infer a specimen calibration, resolved surface geometry, coatings,
//! subsurface transport or polarization retention.

use core::f64::consts::PI;
use fs_math::det;

const C1: f64 = 0.5 - 2.0 / (3.0 * PI);
const C2: f64 = 2.0 / 3.0 - 28.0 / (15.0 * PI);

/// Native reflection proposal version, independent of the unchanged EON BRDF.
pub const ROUGH_DIFFUSE_SAMPLING_SEMANTICS_VERSION: u32 = 2;

/// Validated isotropic rough diffuse surface, evaluated without allocations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RoughDiffuseSurface {
    roughness: f64,
    a: f64,
    b: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RoughDiffuseError {
    /// Roughness must be finite and in [0,1].
    Roughness,
    /// Single-scattering albedo must be finite and in [0,1].
    Albedo,
    /// Cosines or the projected directional dot product are inconsistent.
    Direction,
    /// The requested grazing BRDF exceeds the representable numerical range.
    NumericalRange,
}

impl core::fmt::Display for RoughDiffuseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            Self::Roughness => "rough diffuse roughness must be in [0,1]",
            Self::Albedo => "rough diffuse facet albedo must be in [0,1]",
            Self::Direction => "invalid rough diffuse directions",
            Self::NumericalRange => "rough diffuse BRDF exceeds numerical range",
        })
    }
}

impl core::error::Error for RoughDiffuseError {}

impl RoughDiffuseSurface {
    /// The exact Lambertian limit, retaining the same spectral facet albedo.
    pub const SMOOTH: Self = Self {
        roughness: 0.0,
        a: 1.0,
        b: 0.0,
    };

    pub fn try_new(roughness: f64) -> Result<Self, RoughDiffuseError> {
        if !roughness.is_finite() || !(0.0..=1.0).contains(&roughness) {
            return Err(RoughDiffuseError::Roughness);
        }
        if roughness == 0.0 {
            return Ok(Self::SMOOTH);
        }
        let a = 1.0 / (1.0 + C1 * roughness);
        Ok(Self {
            roughness,
            a,
            b: roughness * a,
        })
    }

    pub const fn roughness(self) -> f64 {
        self.roughness
    }

    // At grazing views the positive FON term has nonzero projected flux as
    // mu_i tends to zero. A cosine-only proposal therefore has large weights.
    // Bound that term by (rho*A/pi)*roughness*sin(theta_o)*max(cos(phi),0).
    // Mix the normalized densities mu_i/pi and max(cos(phi),0)/2 using their
    // envelope masses 1 and 2*roughness*sin(theta_o)/pi. The cosine component
    // retains support for the full BRDF, including its multiple-return term.
    pub(crate) fn grazing_proposal_probability(self, mu_o: f64) -> f64 {
        let sine = ((1.0 - mu_o) * (1.0 + mu_o)).max(0.0).sqrt();
        let mass = (2.0 / PI) * self.roughness * sine;
        mass / (1.0 + mass)
    }

    /// Matching solid-angle density, conditional on the original outgoing
    /// direction. Reciprocity of the BRDF does not make this PDF symmetric.
    pub(crate) fn pdf_admitted(self, mu_o: f64, mu_i: f64, projected_dot: f64) -> f64 {
        if mu_o <= 0.0 || mu_i <= 0.0 { return 0.0; }
        let probability = self.grazing_proposal_probability(mu_o);
        let cosine_pdf = mu_i / PI;
        if probability == 0.0 { return cosine_pdf; }
        let radius = (((1.0 - mu_o) * (1.0 + mu_o))
            * ((1.0 - mu_i) * (1.0 + mu_i))).max(0.0).sqrt();
        // At the isolated pole azimuth has no meaning; choose one common
        // density convention for sampling, connections and reverse queries.
        let azimuth_cosine = if radius > 0.0 {
            (projected_dot / radius).clamp(0.0, 1.0)
        } else { 0.0 };
        (1.0 - probability) * cosine_pdf + probability * (0.5 * azimuth_cosine)
    }

    /// Unit local direction with +X along the outgoing tangential direction.
    /// Admitted samples lie in [0,1); branch selection reuses the first
    /// coordinate, so the native paths retain exactly two random dimensions.
    pub(crate) fn sample_local(self, mu_o: f64, samples: [f64; 2]) -> [f64; 3] {
        let probability = self.grazing_proposal_probability(mu_o);
        let (cosine, radius, sine_phi, cosine_phi) = if samples[0] < probability {
            let cosine = 1.0 - samples[0] / probability;
            let radius = ((1.0 - cosine) * (1.0 + cosine)).sqrt();
            // Uniform elevation and a cosine azimuth on [-pi/2,pi/2].
            // Its inverse CDF needs no trigonometric evaluation.
            let sine_phi = 2.0 * samples[1] - 1.0;
            let cosine_phi = 2.0 * (samples[1] * (1.0 - samples[1])).sqrt();
            (cosine, radius, sine_phi, cosine_phi)
        } else {
            let radial_squared = (samples[0] - probability) / (1.0 - probability);
            let phi = core::f64::consts::TAU * samples[1];
            ((1.0 - radial_squared).sqrt(), radial_squared.sqrt(), det::sin(phi), det::cos(phi))
        };
        [radius * cosine_phi, radius * sine_phi, cosine]
    }

    /// BRDF [1/sr] for unit directions pointing away from the surface.
    /// `projected_dot = wi.wo - (n.wi)(n.wo)` is their tangential dot product,
    /// not their azimuth cosine. Both normal cosines lie in [0,1]. At exactly
    /// zero incoming or outgoing normal flux the measure-zero value is zero.
    pub fn brdf(
        self,
        incident_cosine: f64,
        outgoing_cosine: f64,
        projected_dot: f64,
        facet_albedo: f64,
    ) -> Result<f64, RoughDiffuseError> {
        admit_fraction(facet_albedo, RoughDiffuseError::Albedo)?;
        admit_fraction(incident_cosine, RoughDiffuseError::Direction)?;
        admit_fraction(outgoing_cosine, RoughDiffuseError::Direction)?;
        let bound = ((1.0 - incident_cosine * incident_cosine)
            * (1.0 - outgoing_cosine * outgoing_cosine))
            .max(0.0)
            .sqrt();
        if !projected_dot.is_finite() || projected_dot.abs() > bound + 16.0 * f64::EPSILON {
            return Err(RoughDiffuseError::Direction);
        }
        let value = self.eval_admitted(
            incident_cosine,
            outgoing_cosine,
            projected_dot.clamp(-bound, bound),
            facet_albedo,
        );
        if value.is_finite() {
            Ok(value)
        } else {
            Err(RoughDiffuseError::NumericalRange)
        }
    }

    /// Reflected fraction under collimated illumination, integrated over all
    /// reflected directions. It is exactly one for white facets and zero for
    /// black facets; colored facets lose power at each additional return.
    pub fn directional_reflectance(
        self,
        cosine: f64,
        facet_albedo: f64,
    ) -> Result<f64, RoughDiffuseError> {
        admit_fraction(cosine, RoughDiffuseError::Direction)?;
        admit_fraction(facet_albedo, RoughDiffuseError::Albedo)?;
        if facet_albedo == 0.0 || facet_albedo == 1.0 || self.roughness == 0.0 {
            return Ok(facet_albedo);
        }
        let average = self.single_scatter_average();
        let deficit = self.b * angular_deficit(cosine);
        let denominator = (1.0 - facet_albedo) + facet_albedo * average;
        Ok(facet_albedo * (1.0 - (1.0 - facet_albedo) * deficit / denominator))
    }

    /// Fraction absorbed from a collimated incident direction. By reciprocal
    /// Kirchhoff balance this is also directional emissivity for an opaque
    /// surface in local thermal equilibrium. This law supplies no temperature
    /// or emission spectrum itself.
    ///
    /// Factor the facet absorptance before evaluating the multiple returns:
    /// subtracting a near-unit reflected fraction would erase weak absorption.
    /// At zero cosine this returns the grazing limit, not a finite flux from
    /// an exactly tangent ray.
    pub fn directional_absorptance(
        self,
        cosine: f64,
        facet_albedo: f64,
    ) -> Result<f64, RoughDiffuseError> {
        admit_fraction(cosine, RoughDiffuseError::Direction)?;
        admit_fraction(facet_albedo, RoughDiffuseError::Albedo)?;
        let absorption = 1.0 - facet_albedo;
        if absorption == 0.0 || facet_albedo == 0.0 || self.roughness == 0.0 {
            return Ok(absorption);
        }
        let denominator = absorption + facet_albedo * self.single_scatter_average();
        let deficit = self.b * angular_deficit(cosine);
        Ok(absorption * (1.0 + facet_albedo * deficit / denominator))
    }

    /// Projected-hemisphere mean absorptance, normalized by pi. For the same
    /// opaque reciprocal surface this is the hemispherical emissivity at the
    /// original wavelength, so spectral radiant exitance is pi times this
    /// fraction times the equilibrium blackbody radiance.
    ///
    /// The analytic multiple-return sum is (1-rho)/(1-rho+rho*<E_single>).
    /// It retains weak absorption and has exact smooth, black and white limits.
    pub fn hemispherical_absorptance(self, facet_albedo: f64) -> Result<f64, RoughDiffuseError> {
        admit_fraction(facet_albedo, RoughDiffuseError::Albedo)?;
        let absorption = 1.0 - facet_albedo;
        if absorption == 0.0 || facet_albedo == 0.0 || self.roughness == 0.0 {
            return Ok(absorption);
        }
        Ok(absorption / (absorption + facet_albedo * self.single_scatter_average()))
    }

    pub(crate) fn eval_admitted(self, mu_i: f64, mu_o: f64, projected_dot: f64, rho: f64) -> f64 {
        if mu_i <= 0.0 || mu_o <= 0.0 || rho == 0.0 {
            return 0.0;
        }
        if self.roughness == 0.0 {
            return rho / PI;
        }
        let s = projected_dot.clamp(-1.0, 1.0);
        let angular = if s > 0.0 { s / mu_i.max(mu_o) } else { s };
        let single = (rho / PI) * self.a * (1.0 + self.roughness * angular);
        let average = self.single_scatter_average();
        let multiple_albedo = rho * rho * average / ((1.0 - rho) + rho * average);
        // Factor 1-E(mu)=B*d(mu) and 1-<E>=B*(C1-C2) before division.
        // No epsilon floor or subtraction of nearly unit albedos is needed,
        // including arbitrarily small roughness and the exact smooth limit.
        let multiple = (multiple_albedo / PI)
            * self.b
            * (angular_deficit(mu_i) * angular_deficit(mu_o) / (C1 - C2));
        single + multiple
    }

    fn single_scatter_average(self) -> f64 {
        1.0 - self.b * (C1 - C2)
    }
}

fn admit_fraction(value: f64, error: RoughDiffuseError) -> Result<(), RoughDiffuseError> {
    if value.is_finite() && (0.0..=1.0).contains(&value) {
        Ok(())
    } else {
        Err(error)
    }
}

/// C1-G(mu)/pi from the exact projected FON angular integral (eq. 13).
/// Rewrite acos(mu)=pi/2-asin(mu) and rationalize 1-sin(theta). This keeps
/// the physical O(mu) grazing deficit instead of subtracting two constants.
fn angular_deficit(mu: f64) -> f64 {
    let sin = ((1.0 - mu) * (1.0 + mu)).max(0.0).sqrt();
    let one_minus_sin = mu * mu / (1.0 + sin);
    C1 * one_minus_sin
        + (sin / PI) * (det::asin(mu) - mu * one_minus_sin * (sin + 2.0) / (3.0 * (1.0 + sin)))
}

#[cfg(test)]
#[path = "rough_diffuse/tests.rs"]
mod tests;

//! A transparent smooth dielectric above an ideal diffuse, depolarizing base.
//!
//! The coating is geometrically collapsed and lossless, with constant absolute
//! phase indices. Its exterior index must match the actual surrounding medium.
//! Either index ordering is supported, including external total internal
//! reflection when the exterior has the larger index. There is no finite optical
//! thickness, coherent-film interference, lateral displacement, roughness,
//! bulk absorption or inferred material spectrum. The substrate's reflectance
//! is supplied at the actual wavelength by the existing spectral owner.
//!
//! Entry/exit Fresnel transmission and Snell etendue are retained. Every base
//! reflection depolarizes, so all repeated internal returns form a geometric
//! series. This is the physical nonlinear smooth-plastic model described by
//! Mitsuba, not an additive diffuse/specular blend or color-preserving remap:
//! <https://mitsuba.readthedocs.io/en/latest/src/generated/plugins_bsdfs.html#smooth-plastic-material-plastic>.

use crate::dielectric::{DielectricError, admitted_dielectric_index, fresnel_dielectric};

/// An admitted smooth, lossless coat and its hemispherical escape probability.
///
/// Indices and all probabilities are dimensionless. Construction performs the
/// bounded angular integration once; reflection queries allocate no memory.
/// The retained value is a converged numerical quadrature, not a certified
/// error enclosure or a measured material calibration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmoothDiffuseCoating {
    exterior_index: f64,
    coating_index: f64,
    escape_probability: f64,
}

impl SmoothDiffuseCoating {
    /// Prepare a coat with both absolute phase indices at least one.
    ///
    /// Original dielectric admission also bounds squared indices. An
    /// unrepresentable escape probability or unconverged angular integral
    /// returns `InvalidInterface`; no perfect trapping is manufactured.
    pub fn try_new(exterior_index: f64, coating_index: f64) -> Result<Self, DielectricError> {
        if !admitted_dielectric_index(exterior_index)
            || !admitted_dielectric_index(coating_index) {
            return Err(DielectricError::InvalidInterface);
        }
        let escape_probability = if exterior_index == coating_index { 1.0 } else {
            // Always integrate from the lower index into the higher one,
            // whose complete hemisphere has no critical-angle cusp. For a
            // high-index coat, reciprocity then supplies its internal escape
            // cone's etendue Jacobian. For a low-index coat, this integral is
            // already the base's physical escape probability.
            let low = exterior_index.min(coating_index);
            let high = exterior_index.max(coating_index);
            let mut previous = external_transmission(low, high, 1)?;
            let mut converged = None;
            for level in 1..=8 {
                let next = external_transmission(low, high, 1 << level)?;
                if (next - previous).abs() <= 2.0e-13 * next {
                    converged = Some(next);
                    break;
                }
                previous = next;
            }
            let mean = converged.ok_or(DielectricError::InvalidInterface)?;
            if exterior_index < coating_index {
                let ratio = exterior_index / coating_index;
                (mean * ratio) * ratio
            } else {
                mean
            }
        };
        if !escape_probability.is_finite()
            || escape_probability <= 0.0 || escape_probability > 1.0 {
            return Err(DielectricError::InvalidInterface);
        }
        Ok(Self { exterior_index, coating_index, escape_probability })
    }

    pub const fn exterior_index(self) -> f64 { self.exterior_index }
    pub const fn coating_index(self) -> f64 { self.coating_index }

    /// Probability that unpolarized Lambertian flux at the base escapes the
    /// internal interface on its next encounter, including its escape cone.
    pub const fn escape_probability(self) -> f64 { self.escape_probability }

    pub(crate) const fn is_matched(self) -> bool { self.exterior_index == self.coating_index }

    /// Discrete mirror and diffuse proposal masses, before substrate loss.
    /// Their continuous diffuse density additionally includes the normalized
    /// cosine proposal on the external escape cone.
    pub fn branch_probabilities(self, exterior_cosine: f64) -> Result<[f64; 2], DielectricError> {
        Ok(fresnel_dielectric(exterior_cosine, self.exterior_index, self.coating_index)?
            .event_probabilities())
    }

    /// Continuous unpolarized reflection BRDF in inverse steradians.
    /// The separate specular delta retains the native dielectric reflection.
    pub fn diffuse_brdf(self, incident_cosine: f64, outgoing_cosine: f64, reflectance: f64)
        -> Result<f64, DielectricError> {
        let incoming = fresnel_dielectric(incident_cosine, self.exterior_index, self.coating_index)?;
        let outgoing = fresnel_dielectric(outgoing_cosine, self.exterior_index, self.coating_index)?;
        if incoming.total_internal_reflection || outgoing.total_internal_reflection {
            self.return_denominator(reflectance)?;
            return Ok(0.0);
        }
        Ok((self.diffuse_weight(reflectance)? / core::f64::consts::PI)
            * incoming.transmittance * outgoing.transmittance)
    }

    /// Total reflected fraction for unpolarized collimated illumination.
    /// Includes the mirror and every diffuse internal return.
    pub fn directional_reflectance(self, incident_cosine: f64, reflectance: f64)
        -> Result<f64, DielectricError> {
        let interface = fresnel_dielectric(incident_cosine, self.exterior_index, self.coating_index)?;
        let denominator = self.return_denominator(reflectance)?;
        let diffuse = if reflectance == 1.0 { 1.0 }
            else { reflectance * self.escape_probability / denominator };
        Ok(interface.reflectance + interface.transmittance * diffuse)
    }

    /// The diffuse Mueller coefficient before either interface transmission.
    /// Cosine sampling cancels pi; connection evaluation retains that divisor.
    pub(crate) fn diffuse_weight(self, reflectance: f64) -> Result<f64, DielectricError> {
        let denominator = self.return_denominator(reflectance)?;
        let ratio = self.exterior_index / self.coating_index;
        let result = reflectance * ((ratio * ratio) / denominator);
        if result.is_finite() { Ok(result) } else { Err(DielectricError::InvalidInterface) }
    }

    /// Projected solid-angle fraction occupied by rays escaping the coat.
    /// Multiplying a cosine sampler's squared radius by this value samples
    /// that cone without rejection or a critical-angle subtraction.
    pub(crate) fn diffuse_cone_fraction(self) -> f64 {
        if self.coating_index >= self.exterior_index { 1.0 } else {
            let ratio = self.coating_index / self.exterior_index;
            ratio * ratio
        }
    }

    /// Conditional directional density; a low-index coating cannot return
    /// diffuse light outside its Snell escape cone in the exterior medium.
    pub(crate) fn diffuse_direction_pdf(self, cosine: f64) -> Result<f64, DielectricError> {
        if !cosine.is_finite() || !(0.0..=1.0).contains(&cosine) {
            return Err(DielectricError::InvalidInterface);
        }
        if self.coating_index < self.exterior_index
            && fresnel_dielectric(cosine, self.exterior_index, self.coating_index)?
                .total_internal_reflection {
            return Ok(0.0);
        }
        Ok((cosine / core::f64::consts::PI) / self.diffuse_cone_fraction())
    }

    /// Cancel the cone proposal's etendue analytically before multiplying
    /// powers. Forming the BRDF and dividing by its narrow-cone density first
    /// would introduce large, mutually cancelling intermediate coefficients.
    pub(crate) fn sampled_diffuse_weight(self, reflectance: f64) -> Result<f64, DielectricError> {
        if self.coating_index >= self.exterior_index {
            self.diffuse_weight(reflectance)
        } else {
            Ok(reflectance / self.return_denominator(reflectance)?)
        }
    }

    fn return_denominator(self, reflectance: f64) -> Result<f64, DielectricError> {
        if !(0.0..=1.0).contains(&reflectance) { return Err(DielectricError::InvalidAbsorption); }
        // 1-rho*(1-escape), written without cancellation near a white base.
        Ok((1.0 - reflectance) + reflectance * self.escape_probability)
    }
}

// Positive abscissae and weights of the standard eight-point Gauss-Legendre
// rule. Composite refinement is confined to construction, never a ray event.
const GAUSS_8: [(f64, f64); 4] = [
    (0.1834346424956498, 0.3626837833783620),
    (0.5255324099163290, 0.3137066458778873),
    (0.7966664774136267, 0.2223810344533745),
    (0.9602898564975363, 0.1012285362903763),
];

fn external_transmission(exterior: f64, coating: f64, panels: usize)
    -> Result<f64, DielectricError> {
    let half_width = 0.5 / panels as f64;
    let mut sum = 0.0;
    let mut correction = 0.0;
    for panel in 0..panels {
        let midpoint = (panel as f64 + 0.5) / panels as f64;
        for (abscissa, weight) in GAUSS_8 {
            for sign in [-1.0, 1.0] {
                let cosine = midpoint + sign * half_width * abscissa;
                let transmitted = fresnel_dielectric(cosine, exterior, coating)?.transmittance;
                let value = weight * cosine * transmitted - correction;
                let next = sum + value;
                correction = (next - sum) - value;
                sum = next;
            }
        }
    }
    Ok(2.0 * half_width * sum)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Independent textbook s/p fields, with Snell geometry expressed directly.
    fn transmission(cosine: f64, eta: f64) -> f64 {
        if eta == 1.0 { return 1.0; }
        let transmitted_squared = 1.0 - (1.0 - cosine * cosine) / (eta * eta);
        if transmitted_squared <= 0.0 { return 0.0; }
        let transmitted = transmitted_squared.sqrt();
        let s = (cosine - eta * transmitted) / (cosine + eta * transmitted);
        let p = (eta * cosine - transmitted) / (eta * cosine + transmitted);
        1.0 - 0.5 * (s * s + p * p)
    }

    // Composite Simpson uses neither the constructor's nodes nor native
    // Fresnel code. Its larger fixed grid is only an independent test oracle.
    fn integrate(mut f: impl FnMut(f64) -> f64) -> f64 {
        const N: usize = 65536;
        let mut sum = f(0.0) + f(1.0);
        for i in 1..N { sum += if i % 2 == 0 { 2.0 } else { 4.0 } * f(i as f64 / N as f64); }
        sum / (3.0 * N as f64)
    }

    #[test]
    fn g1_furnace_and_internal_returns_match_independent_flux_integrals() {
        for eta in [1.0, 1.01, 1.3, 1.5, 2.4] {
            let coating = SmoothDiffuseCoating::try_new(1.0, eta).unwrap();
            let external_escape = integrate(|mu| 2.0 * mu * transmission(mu, eta));
            assert!((coating.escape_probability - external_escape / (eta * eta)).abs() < 2e-11);
            for cosine in [0.03, 0.3, 0.8, 1.0] {
                let entry = transmission(cosine, eta);
                for rho in [0.0, 0.2, 0.8, 1.0] {
                    let reflected = 1.0 - entry;
                    let scale = coating.diffuse_weight(rho).unwrap();
                    let furnace = reflected + entry * scale * external_escape;
                    let escape = external_escape / (eta * eta);
                    let expected = reflected + entry * rho * escape / ((1.0-rho) + rho*escape);
                    assert!((furnace - expected).abs() < 2e-11);
                    assert!((coating.directional_reflectance(cosine, rho).unwrap() - expected).abs() < 2e-11);
                    assert!(furnace >= 0.0 && furnace <= 1.0 + 2e-11);
                    if rho == 1.0 { assert!((furnace - 1.0).abs() < 2e-11); }
                }
            }
        }
    }

    #[test]
    fn g3_reciprocity_matched_indices_and_absolute_index_scaling() {
        let coating = SmoothDiffuseCoating::try_new(1.0, 1.5).unwrap();
        let scaled = SmoothDiffuseCoating::try_new(1.4, 2.1).unwrap();
        for incoming in [0.01, 0.2, 0.7, 1.0] {
            for outgoing in [0.03, 0.4, 0.9, 1.0] {
                let value = coating.diffuse_brdf(incoming, outgoing, 0.6).unwrap();
                let reverse = coating.diffuse_brdf(outgoing, incoming, 0.6).unwrap();
                assert!((value - reverse).abs() < 4.0 * f64::EPSILON);
                assert!((value - scaled.diffuse_brdf(incoming, outgoing, 0.6).unwrap()).abs() < 2e-13);
                let matched = SmoothDiffuseCoating::try_new(1.33, 1.33).unwrap();
                assert_eq!(matched.diffuse_brdf(incoming, outgoing, 0.6).unwrap(), 0.6 / core::f64::consts::PI);
                assert_eq!(matched.branch_probabilities(incoming).unwrap(), [0.0, 1.0]);
            }
        }
    }

    #[test]
    fn g1_low_index_coats_retain_escape_cone_and_all_internal_returns() {
        for (outside, inside) in [(1.5, 1.0), (2.4, 1.3), (4.0, 1.01)] {
            let coating = SmoothDiffuseCoating::try_new(outside, inside).unwrap();
            let ratio = inside / outside;
            let cone = ratio * ratio;
            let escape = integrate(|mu| 2.0 * mu * transmission(mu, 1.0 / ratio));
            assert!((coating.escape_probability() - escape).abs() < 2e-11);
            assert_eq!(coating.diffuse_cone_fraction(), cone);
            for sine_squared in [0.0, 0.2, 0.9, 1.1] {
                let incoming = (1.0 - cone * sine_squared).sqrt();
                let entry = transmission(incoming, ratio);
                for rho in [0.0, 0.3, 0.9, 1.0] {
                    // Integrate over the ORIGINAL internal hemisphere. Its
                    // projected etendue maps to the external escape cone;
                    // this avoids reusing the production sampling measure.
                    let diffuse = integrate(|internal_mu| {
                        let external_mu = (1.0 - cone * (1.0 - internal_mu*internal_mu)).sqrt();
                        2.0 * core::f64::consts::PI * cone * internal_mu
                            * coating.diffuse_brdf(incoming, external_mu, rho).unwrap()
                    });
                    let expected = (1.0-entry) + entry*rho*escape/((1.0-rho)+rho*escape);
                    assert!(((1.0-entry)+diffuse - expected).abs() < 2e-10);
                    assert!((coating.directional_reflectance(incoming, rho).unwrap()-expected).abs() < 2e-11);
                    if rho == 1.0 { assert!(((1.0-entry)+diffuse - 1.0).abs() < 2e-10); }
                }
            }
        }
    }

    #[test]
    fn g3_immersed_coat_reciprocity_and_total_internal_reflection() {
        let coating = SmoothDiffuseCoating::try_new(1.5, 1.0).unwrap();
        let scaled = SmoothDiffuseCoating::try_new(3.0, 2.0).unwrap();
        for incident in [0.2, 0.8, 0.9, 1.0] {
            for outgoing in [0.3, 0.78, 0.95, 1.0] {
                let forward = coating.diffuse_brdf(incident, outgoing, 0.6).unwrap();
                let reverse = coating.diffuse_brdf(outgoing, incident, 0.6).unwrap();
                assert!((forward-reverse).abs() < 4e-16);
                assert_eq!(forward, scaled.diffuse_brdf(incident, outgoing, 0.6).unwrap());
            }
        }
        assert_eq!(coating.branch_probabilities(0.5).unwrap(), [1.0, 0.0]);
        assert_eq!(coating.directional_reflectance(0.5, 0.0).unwrap(), 1.0);
        assert_eq!(coating.diffuse_brdf(0.5, 1.0, 1.0).unwrap(), 0.0);
        assert_eq!(coating.diffuse_direction_pdf(0.5).unwrap(), 0.0);
        assert_eq!(coating.diffuse_direction_pdf(1.0).unwrap(),
            (1.0 / core::f64::consts::PI) / (1.0_f64 / 1.5).powi(2));
    }

    #[test]
    fn g0_coating_domains_and_black_substrate_remain_explicit() {
        for (outside, inside) in [(0.0, 1.5), (1.5, 0.9), (f64::NAN, 1.5), (1.0, f64::INFINITY)] {
            assert!(SmoothDiffuseCoating::try_new(outside, inside).is_err());
        }
        let coating = SmoothDiffuseCoating::try_new(1.0, 1.5).unwrap();
        assert_eq!(coating.diffuse_brdf(0.5, 0.7, 0.0).unwrap(), 0.0);
        assert!((coating.directional_reflectance(1.0, 0.0).unwrap() - 0.04).abs() < 1e-15);
        for rho in [-0.1, 1.01, f64::NAN] { assert!(coating.diffuse_brdf(0.5, 0.7, rho).is_err()); }
        assert!(coating.diffuse_brdf(-0.1, 0.7, 0.5).is_err());
    }
}

//! A transparent smooth dielectric above an ideal diffuse, depolarizing base.
//!
//! The coating is geometrically collapsed and lossless, with constant absolute
//! phase indices. Its exterior index must match the actual surrounding medium.
//! The interior/exterior ratio is at least one. There is no finite optical
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
    /// Prepare a coat with `1 <= exterior_index <= coating_index`.
    ///
    /// Original dielectric admission also bounds squared indices. An
    /// unrepresentable escape probability or unconverged angular integral
    /// returns `InvalidInterface`; no perfect trapping is manufactured.
    pub fn try_new(exterior_index: f64, coating_index: f64) -> Result<Self, DielectricError> {
        if !admitted_dielectric_index(exterior_index)
            || !admitted_dielectric_index(coating_index)
            || coating_index < exterior_index {
            return Err(DielectricError::InvalidInterface);
        }
        let escape_probability = if exterior_index == coating_index { 1.0 } else {
            // Reciprocity maps the internal escape cone to the complete
            // external hemisphere. Integrating there avoids the internal TIR
            // cusp. The two index ratios are the Snell etendue Jacobian.
            let mut previous = external_transmission(exterior_index, coating_index, 1)?;
            let mut converged = None;
            for level in 1..=8 {
                let next = external_transmission(exterior_index, coating_index, 1 << level)?;
                if (next - previous).abs() <= 2.0e-13 * next {
                    converged = Some(next);
                    break;
                }
                previous = next;
            }
            let mean = converged.ok_or(DielectricError::InvalidInterface)?;
            let ratio = exterior_index / coating_index;
            (mean * ratio) * ratio
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
    /// Their continuous diffuse density is `mass[1] * cos(theta_out) / pi`.
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
        let transmitted = (1.0 - (1.0 - cosine * cosine) / (eta * eta)).sqrt();
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
    fn g0_coating_domains_and_black_substrate_remain_explicit() {
        for (outside, inside) in [(0.0, 1.5), (1.5, 1.0), (f64::NAN, 1.5), (1.0, f64::INFINITY)] {
            assert!(SmoothDiffuseCoating::try_new(outside, inside).is_err());
        }
        let coating = SmoothDiffuseCoating::try_new(1.0, 1.5).unwrap();
        assert_eq!(coating.diffuse_brdf(0.5, 0.7, 0.0).unwrap(), 0.0);
        assert!((coating.directional_reflectance(1.0, 0.0).unwrap() - 0.04).abs() < 1e-15);
        for rho in [-0.1, 1.01, f64::NAN] { assert!(coating.diffuse_brdf(0.5, 0.7, rho).is_err()); }
        assert!(coating.diffuse_brdf(-0.1, 0.7, 0.5).is_err());
    }
}

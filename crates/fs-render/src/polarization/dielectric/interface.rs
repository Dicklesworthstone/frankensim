//! Polarized smooth-interface transport using the renderer's existing glass,
//! Fresnel admission, refraction geometry and radiance sampling conventions.
//!
//! Both Stokes frames have +Q along s (perpendicular to the incidence plane),
//! with p chosen so (s,p,propagation direction) is right-handed. Consequently
//! reflection reverses the p axis at normal incidence. We use
//! U=2 Re(E_s conj(E_p)), V=2 Im(E_s conj(E_p)), and exp(-i omega t).
//! Rotate into/out of these frames at every vertex; unrelated ray frames must
//! never be composed directly. The returned transmission matrix transfers
//! *power*, not radiance; the sampling adapter includes the legacy eta² factor.
//!
//! Only lossless isotropic interfaces are admitted. Existing Beer-Lambert
//! absorption is a separate bulk effect, not a complex interface index.
//! Formula reference: Mitsuba's polarization/Fresnel conventions:
//! <https://mitsuba.readthedocs.io/en/stable/src/key_topics/polarization.html>.

use fs_geom::Vec3;

use crate::dielectric::{
    BeerLambertAbsorption, DielectricError, DielectricEvent, DielectricIor, FresnelDielectric,
    fresnel_dielectric, sample_smooth_dielectric,
};
use crate::polarization::{Mueller, Stokes};

/// Power-normalized Fresnel Mueller matrices in the interface s/p frames.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarizedFresnel {
    /// Reflection, including relative phase under total internal reflection.
    pub reflection: Mueller,
    /// Transmission of ray power; zero for total internal reflection.
    pub transmission: Mueller,
    /// Unpolarized reflection and the native Snell geometry/boundary policy.
    pub scalar: FresnelDielectric,
}

pub(super) fn diagonal(s_power: f64, p_power: f64, cross_re: f64, cross_im: f64) -> Mueller {
    let mean = 0.5 * s_power + 0.5 * p_power;
    let difference = 0.5 * s_power - 0.5 * p_power;
    Mueller([
        [mean, difference, 0.0, 0.0],
        [difference, mean, 0.0, 0.0],
        [0.0, 0.0, cross_re, -cross_im],
        [0.0, 0.0, cross_im, cross_re],
    ])
}

/// (a-i b)/(a+i b), evaluated without complex allocation or atan/trig.
fn tir_amplitude(a: f64, b: f64) -> (f64, f64) {
    let denominator = a * a + b * b;
    ((a * a - b * b) / denominator, -2.0 * a * b / denominator)
}

/// Exact s/p dielectric response for the existing `[1,3]` index domain.
///
/// This deliberately calls the scalar admission and Snell calculation, so
/// equal-index/grazing/critical-angle policy cannot diverge from the tracer.
/// No allocation or transcendental phase evaluation occurs per interface.
pub fn fresnel_dielectric_polarized(
    incident_cosine: f64,
    eta_i: f64,
    eta_t: f64,
) -> Result<PolarizedFresnel, DielectricError> {
    let scalar = fresnel_dielectric(incident_cosine, eta_i, eta_t)?;
    let zero = Mueller([[0.0; 4]; 4]);
    if eta_i == eta_t {
        return Ok(PolarizedFresnel {
            reflection: zero,
            transmission: Mueller::IDENTITY,
            scalar,
        });
    }
    if scalar.total_internal_reflection {
        let ratio = eta_i / eta_t;
        let gamma = (ratio * ratio * (1.0 - incident_cosine * incident_cosine)
            - 1.0).max(0.0).sqrt();
        let (sr, si) = tir_amplitude(eta_i * incident_cosine, eta_t * gamma);
        let (pr, pi) = tir_amplitude(eta_t * incident_cosine, eta_i * gamma);
        // r_s conj(r_p): retain the phase difference even though both powers
        // are one. Replacing this by identity loses elliptical polarization.
        return Ok(PolarizedFresnel {
            reflection: diagonal(1.0, 1.0, sr * pr + si * pi, si * pr - sr * pi),
            transmission: zero,
            scalar,
        });
    }
    let ci = incident_cosine;
    let ct = scalar.transmitted_cosine;
    let ds = eta_i * ci + eta_t * ct;
    let dp = eta_t * ci + eta_i * ct;
    let rs = (eta_i * ci - eta_t * ct) / ds;
    let rp = (eta_t * ci - eta_i * ct) / dp;
    // Power-normalized transmission, not squared E-field transmission.
    // The direct expression avoids cancellation in 1-r² near grazing.
    let numerator = 4.0 * eta_i * eta_t * ci * ct;
    let ts = (numerator / (ds * ds)).clamp(0.0, 1.0);
    let tp = (numerator / (dp * dp)).clamp(0.0, 1.0);
    Ok(PolarizedFresnel {
        reflection: diagonal(rs * rs, rp * rp, rs * rp, 0.0),
        transmission: diagonal(ts, tp, ts.sqrt() * tp.sqrt(), 0.0),
        scalar,
    })
}

/// Evaluate the same phase-index laws already used by the spectral tracer.
/// Indices are in incident/transmitted order; wavelength is vacuum nanometres.
pub fn spectral_fresnel_dielectric_polarized(
    incident_cosine: f64,
    incident: impl Into<DielectricIor>,
    transmitted: impl Into<DielectricIor>,
    wavelength_nm: f64,
) -> Result<PolarizedFresnel, DielectricError> {
    spectral_fresnel_dielectric_polarized_borrowed(
        incident_cosine,
        &incident.into(),
        &transmitted.into(),
        wavelength_nm,
    )
}

// A spectral filter converts each input law once, then borrows it across bins.
// This keeps inline index tables out of the per-bin conversion path.
pub(in crate::polarization) fn spectral_fresnel_dielectric_polarized_borrowed(
    incident_cosine: f64,
    incident: &DielectricIor,
    transmitted: &DielectricIor,
    wavelength_nm: f64,
) -> Result<PolarizedFresnel, DielectricError> {
    fresnel_dielectric_polarized(
        incident_cosine,
        incident.eval(wavelength_nm)?,
        transmitted.eval(wavelength_nm)?,
    )
}

/// A polarized counterpart of the existing smooth-dielectric Monte Carlo sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarizedDielectricSample {
    /// Unit outgoing direction, using the existing away-from-surface convention.
    pub direction: Vec3,
    /// Sampled reflection/transmission event.
    pub event: DielectricEvent,
    /// Scalar Fresnel proposal probability, not polarization-dependent power.
    pub probability: f64,
    /// Mueller throughput divided by proposal probability. Transmission also
    /// includes the existing `(eta_i/eta_t)^2` radiance-mode factor.
    pub radiance_weight: Mueller,
}

/// Reuse the native smooth sampler rather than fork its geometry or sampling.
///
/// The unpolarized proposal has support for both polarization components;
/// dividing the full matrix by that proposal is unbiased for arbitrary input
/// Stokes states. A scalar path remains unchanged and pays no Mueller cost.
pub fn sample_smooth_dielectric_polarized(
    normal: Vec3,
    wo: Vec3,
    eta_i: f64,
    eta_t: f64,
    event_sample: f64,
) -> Result<PolarizedDielectricSample, DielectricError> {
    let sample = sample_smooth_dielectric(normal, wo, eta_i, eta_t, event_sample)?;
    let n = normal.scale(1.0 / normal.dot(normal).sqrt());
    let w = wo.scale(1.0 / wo.dot(wo).sqrt());
    let fresnel = fresnel_dielectric_polarized(n.dot(w).clamp(0.0, 1.0), eta_i, eta_t)?;
    let (matrix, scale) = match sample.event {
        DielectricEvent::Reflection => (fresnel.reflection, 1.0 / sample.probability),
        DielectricEvent::Transmission => (
            fresnel.transmission, (eta_i / eta_t).powi(2) / sample.probability,
        ),
    };
    Ok(PolarizedDielectricSample {
        direction: sample.direction,
        event: sample.event,
        probability: sample.probability,
        radiance_weight: Mueller(matrix.0.map(|row| row.map(|value| value * scale))),
    })
}

/// Isotropic bulk absorption attenuates all four Stokes components equally.
/// Uses the native material law; it does not depolarize surviving light.
pub fn attenuate_stokes(
    input: Stokes,
    absorption: BeerLambertAbsorption,
    wavelength_nm: f64,
    distance_m: f64,
) -> Result<Stokes, DielectricError> {
    if Stokes::try_new(input.i, input.q, input.u, input.v).is_none() {
        return Err(DielectricError::InvalidAbsorption);
    }
    let t = absorption.transmittance(wavelength_nm, distance_m)?;
    Ok(Stokes { i: input.i * t, q: input.q * t, u: input.u * t, v: input.v * t })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dielectric::{CauchyIor, DielectricIorSample, TabulatedIor};

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() < 2e-12, "{a} != {b}");
    }

    #[test]
    fn normal_incidence_preserves_reflection_basis_phase() {
        let f = fresnel_dielectric_polarized(1.0, 1.0, 1.5).unwrap();
        close(f.reflection.0[0][0], 0.04);
        close(f.transmission.0[0][0], 0.96);
        close(f.reflection.0[2][2], -0.04);
        close(f.reflection.0[3][3], -0.04);
    }

    #[test]
    fn brewster_reflection_is_s_polarized() {
        let ci = 1.0 / 1.5_f64.hypot(1.0);
        let f = fresnel_dielectric_polarized(ci, 1.0, 1.5).unwrap();
        let unpolarized = f.reflection.apply(Stokes::unpolarized(1.0));
        close(unpolarized.q, unpolarized.i);
        close(f.reflection.apply(Stokes::try_new(1.0, -1.0, 0.0, 0.0).unwrap()).i, 0.0);
    }

    #[test]
    fn tir_retains_analytic_relative_phase() {
        let f = fresnel_dielectric_polarized(0.5_f64.sqrt(), 1.5, 1.0).unwrap();
        assert!(f.scalar.total_internal_reflection);
        let out = f.reflection.apply(Stokes::try_new(1.0, 0.0, 1.0, 0.0).unwrap());
        close(out.i, 1.0);
        close(out.u, 0.8);
        close(out.v, 0.6);
        close(out.degree(), 1.0);
        assert_eq!(f.transmission, Mueller([[0.0; 4]; 4]));
    }

    #[test]
    fn power_closure_physical_states_and_scalar_agreement() {
        for eta_i in [1.0, 1.33, 1.5, 3.0] {
            for eta_t in [1.0, 1.33, 1.5, 3.0] {
                for j in 0..=100 {
                    let ci = f64::from(j) / 100.0;
                    let f = fresnel_dielectric_polarized(ci, eta_i, eta_t).unwrap();
                    close(f.reflection.0[0][0], f.scalar.reflectance);
                    for (q, u, v) in [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
                        (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)] {
                        let input = Stokes::try_new(1.0, q, u, v).unwrap();
                        let r = f.reflection.apply(input);
                        let t = f.transmission.apply(input);
                        close(r.i + t.i, 1.0);
                        assert!(Stokes::try_new(r.i, r.q, r.u, r.v).is_some());
                        assert!(Stokes::try_new(t.i, t.q, t.u, t.v).is_some());
                    }
                }
            }
        }
    }

    #[test]
    fn equal_indices_are_identity_even_at_grazing() {
        for ci in [0.0, 0.5, 1.0] {
            let f = fresnel_dielectric_polarized(ci, 1.5, 1.5).unwrap();
            assert_eq!(f.reflection, Mueller([[0.0; 4]; 4]));
            assert_eq!(f.transmission, Mueller::IDENTITY);
        }
    }

    #[test]
    fn sampling_reuses_native_directions_probabilities_and_radiance_scaling() {
        let n = Vec3::new(0.0, 0.0, 1.0);
        let w = Vec3::new(0.6, 0.0, 0.8);
        for u in [0.0, 0.999] {
            let scalar = sample_smooth_dielectric(n, w, 1.0, 1.5, u).unwrap();
            let polarized = sample_smooth_dielectric_polarized(n, w, 1.0, 1.5, u).unwrap();
            assert_eq!(polarized.direction, scalar.direction);
            assert_eq!(polarized.event, scalar.event);
            assert_eq!(polarized.probability, scalar.probability);
            close(polarized.radiance_weight.apply(Stokes::unpolarized(1.0)).i,
                scalar.radiance_weight);
        }
    }

    #[test]
    fn spectral_material_and_absorption_are_shared_not_reimplemented() {
        let air = CauchyIor::try_constant(1.0).unwrap();
        let glass = crate::dielectric::DielectricGlass::representative_crown();
        for wavelength in [400.0, 550.0, 700.0] {
            let f = spectral_fresnel_dielectric_polarized(0.6, air, glass.ior(), wavelength).unwrap();
            assert_eq!(f.scalar, fresnel_dielectric(0.6, 1.0, glass.ior().eval(wavelength).unwrap()).unwrap());
            let input = Stokes::try_new(3.0, 1.0, -1.0, 1.0).unwrap();
            let out = attenuate_stokes(input, glass.absorption(), wavelength, 0.02).unwrap();
            close(out.i / input.i, glass.absorption().transmittance(wavelength, 0.02).unwrap());
            close(out.degree(), input.degree());
        }
        assert!(fresnel_dielectric_polarized(f64::NAN, 1.0, 1.5).is_err());
        assert!(fresnel_dielectric_polarized(0.5, 0.0, 1.5).is_err());
        assert!(spectral_fresnel_dielectric_polarized(0.5, air, air, 100.0).is_err());
    }

    #[test]
    fn tabulated_phase_indices_preserve_polarized_spectral_power() {
        let air = CauchyIor::try_constant(1.0).unwrap();
        let indices = [1.61, 1.60, 1.59, 1.58, 1.57, 1.575, 1.56, 1.55, 1.54];
        let sampled = TabulatedIor::try_new(core::array::from_fn(|bin| {
            DielectricIorSample::try_new(380.0 + 50.0 * bin as f64, indices[bin]).unwrap()
        }))
        .unwrap();
        let law = DielectricIor::from(sampled);
        let polarized = Stokes::try_new(1.0, 0.3, 0.4, 0.5).unwrap();
        for (bin, index) in indices.into_iter().enumerate() {
            let wavelength = 380.0 + 50.0 * bin as f64;
            let actual = spectral_fresnel_dielectric_polarized(0.6, air, &law, wavelength).unwrap();
            let expected = fresnel_dielectric_polarized(0.6, 1.0, index).unwrap();
            assert_eq!(actual, expected);
            let reflected = actual.reflection.apply(polarized);
            let transmitted = actual.transmission.apply(polarized);
            close(reflected.i + transmitted.i, polarized.i);
        }
        let midpoint = spectral_fresnel_dielectric_polarized(0.6, air, sampled, 405.0).unwrap();
        assert_eq!(
            midpoint,
            fresnel_dielectric_polarized(0.6, 1.0, 0.5 * (indices[0] + indices[1])).unwrap()
        );
        assert!(spectral_fresnel_dielectric_polarized(0.6, air, &law, 379.0).is_err());
    }
}

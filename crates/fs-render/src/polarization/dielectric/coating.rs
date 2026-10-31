//! Coherent transparent films on native dielectric ray boundaries.
//!
//! The original coating owns its visible phase-index law, normal thickness and
//! source assertion. Both surrounding media retain the renderer's real-index
//! interface plus separate Beer-Lambert bulk absorption model. Fields use
//! exp(-i omega t); +Q is s, and reflected p reverses at normal incidence.
//!
//! This local smooth-film model collapses the two parallel film faces onto the
//! geometric boundary. It retains interference and relative s/p phase, not a
//! nanometre ray displacement, group delay or coherence between scene paths.
//! The caller owns the film's coherence and locally planar approximation.
//!
//! Transfer/Poynting conventions: S. J. Byrnes, "Multilayer optical
//! calculations", sections 3, 4.1 and 4.4, <https://arxiv.org/abs/1603.02720>.

use crate::conductor::CoherentDielectricCoating;
use crate::dielectric::{
    DielectricCoatingStack, DielectricError, DielectricEvent, fresnel_dielectric,
    normalize_admitted, reflect, refract, validate_direction_frame,
};
use crate::polarization::Mueller;
use fs_geom::Vec3;

use super::{
    PolarizedDielectricSample, PolarizedFresnel, fresnel_dielectric_polarized,
    sample_smooth_dielectric_polarized,
};

/// Coherent s/p reflection and power-normalized transmission through one
/// original transparent coating. Wavelength is vacuum nm; indices retain the
/// native [1,3] domain. The same film is used when incident/exit media reverse.
///
/// A zero-thickness or index-matched film preserves the entire bare Mueller
/// response. Equal surrounding indices do not remove a contrasting film.
/// Exact grazing for a nontrivial film refuses: incident normal power is zero
/// and the resonant grazing limit is not uniquely specified by that angle.
pub fn fresnel_coated_dielectric_polarized(
    incident_cosine: f64,
    eta_i: f64,
    eta_t: f64,
    wavelength_nm: f64,
    coating: &CoherentDielectricCoating,
) -> Result<PolarizedFresnel, DielectricError> {
    fresnel_dielectric(incident_cosine, eta_i, eta_t)?;
    let film_eta = coating.ior().eval(wavelength_nm)?;
    if coating.thickness_m() == 0.0 || film_eta == eta_i || film_eta == eta_t {
        return fresnel_dielectric_polarized(incident_cosine, eta_i, eta_t);
    }
    // Borrow the original spectral film without fitting or copying its law.
    // The shared solver owns both explicit scalar R/T and signed Mueller flux,
    // including the same native exterior critical-angle classification.
    DielectricCoatingStack::try_new(core::slice::from_ref(coating))?.fresnel_polarized(
        wavelength_nm,
        false,
        incident_cosine,
        eta_i,
        eta_t,
    )
}

/// Smooth coated-dielectric radiance sample using the original direction
/// admission, reflection and Snell refraction geometry. The two positive
/// Mueller powers determine the proposal; the smaller branch occupies the
/// beginning of [0,1), retaining weak reflection AND weak transmission.
///
/// Matrices are divided by the selected proposal. Transmission also includes
/// the existing (eta_i/eta_t)^2 radiance factor. Forward beam transport removes
/// that Jacobian exactly once. A zero-thickness film keeps the bare sampler.
pub fn sample_smooth_coated_dielectric_polarized(
    normal: Vec3,
    wo: Vec3,
    eta_i: f64,
    eta_t: f64,
    wavelength_nm: f64,
    coating: &CoherentDielectricCoating,
    event_sample: f64,
) -> Result<PolarizedDielectricSample, DielectricError> {
    if coating.thickness_m() == 0.0 {
        coating.ior().eval(wavelength_nm)?;
        return sample_smooth_dielectric_polarized(normal, wo, eta_i, eta_t, event_sample);
    }
    validate_direction_frame(normal, wo)?;
    if !event_sample.is_finite() || !(0.0..1.0).contains(&event_sample) {
        return Err(DielectricError::InvalidInterface);
    }
    let normal = normalize_admitted(normal);
    let wo = normalize_admitted(wo);
    let response = fresnel_coated_dielectric_polarized(
        normal.dot(wo).clamp(0.0, 1.0),
        eta_i,
        eta_t,
        wavelength_nm,
        coating,
    )?;
    let reflected = response.reflection.0[0][0];
    let transmitted = response.transmission.0[0][0];
    let total = reflected + transmitted;
    let reflection = if reflected <= transmitted {
        event_sample < reflected / total
    } else {
        event_sample >= transmitted / total
    };
    let (direction, event, probability, matrix, radiance) = if reflection {
        (
            reflect(wo, normal),
            DielectricEvent::Reflection,
            reflected / total,
            response.reflection,
            1.0,
        )
    } else {
        (
            refract(wo, normal, eta_i, eta_t, response.scalar.transmitted_cosine)?,
            DielectricEvent::Transmission,
            transmitted / total,
            response.transmission,
            (eta_i / eta_t).powi(2),
        )
    };
    if !probability.is_finite() || probability <= 0.0 {
        return Err(DielectricError::InvalidInterface);
    }
    // Divide first: a subnormal positive proposal must not overflow a
    // reciprocal that would later multiply an equally small optical matrix.
    let radiance_weight = Mueller(
        matrix
            .0
            .map(|row| row.map(|value| (value / probability) * radiance)),
    );
    Ok(PolarizedDielectricSample {
        direction,
        event,
        probability,
        radiance_weight,
    })
}

#[cfg(test)]
mod tests;

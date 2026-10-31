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

use fs_geom::Vec3;
use fs_math::{c64::C64, det};

use crate::conductor::{CoherentDielectricCoating, coating_propagation};
use crate::dielectric::{
    DielectricError, DielectricEvent, fresnel_dielectric, normalize_admitted, reflect, refract,
    validate_direction_frame,
};
use crate::polarization::Mueller;

use super::interface::diagonal;
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
    let mut scalar = fresnel_dielectric(incident_cosine, eta_i, eta_t)?;
    let film_eta = coating.ior().eval(wavelength_nm)?;
    if coating.thickness_m() == 0.0 || film_eta == eta_i || film_eta == eta_t {
        return fresnel_dielectric_polarized(incident_cosine, eta_i, eta_t);
    }
    if incident_cosine <= 0.0 {
        return Err(DielectricError::InvalidInterface);
    }
    let q0 = eta_i * incident_cosine;
    let film_squared = film_eta * film_eta;
    let q1_squared = normal_squared(film_eta, eta_i, incident_cosine);
    let q2 = if scalar.total_internal_reflection {
        C64::new(
            0.0,
            det::sqrt((-normal_squared(eta_t, eta_i, incident_cosine)).max(0.0)),
        )
    } else {
        C64::from_re(eta_t * scalar.transmitted_cosine)
    };
    let phase_scale = core::f64::consts::TAU * (coating.thickness_m() * 1e9) / wavelength_nm;
    // This is the SAME critical-angle-safe bounded propagation used by the
    // conductor's transparent film. Its evanescent common scale must also be
    // carried into transmission; it cancels only from reflection.
    let (cosine, sine_over_q) = coating_propagation(q1_squared, phase_scale);
    let attenuation = if q1_squared < 0.0 {
        det::exp(-phase_scale * det::sqrt(-q1_squared))
    } else {
        1.0
    };
    let propagating = !scalar.total_internal_reflection;
    let (rs, ts) = amplitudes(
        q0,
        q2,
        cosine,
        sine_over_q,
        q1_squared * sine_over_q,
        attenuation,
        propagating,
    )?;
    // H-dual admittance gives the existing rotating-p basis directly and
    // stays finite when the exit's normal wave number is exactly zero.
    let (rp, tp) = if incident_cosine == 1.0 {
        (-rs, ts)
    } else {
        amplitudes(
            incident_cosine / eta_i,
            q2.scale(1.0 / (eta_t * eta_t)),
            cosine,
            film_squared * sine_over_q,
            q1_squared * sine_over_q / film_squared,
            attenuation,
            propagating,
        )?
    };
    let reflection = matrix(rs, rp);
    let transmission = matrix(ts, tp);
    scalar.reflectance = reflection.0[0][0].clamp(0.0, 1.0);
    Ok(PolarizedFresnel {
        reflection,
        transmission,
        scalar,
    })
}

// Preserve grazing c^2 without subtracting two rounded equal squares.
fn normal_squared(index: f64, incident: f64, cosine: f64) -> f64 {
    if cosine > 0.5 {
        index * index - incident * incident * (1.0 - cosine) * (1.0 + cosine)
    } else {
        (index - incident) * (index + incident) + (incident * cosine).powi(2)
    }
}

fn amplitudes(
    incident: f64,
    exit: C64,
    cosine: f64,
    upper: f64,
    lower: f64,
    attenuation: f64,
    propagating: bool,
) -> Result<(C64, C64), DielectricError> {
    let field = C64::from_re(cosine) - C64::new(0.0, 1.0) * exit.scale(upper);
    let flux = exit.scale(cosine) - C64::new(0.0, lower);
    let denominator = field.scale(incident) + flux;
    let reflected = (field.scale(incident) - flux) / denominator;
    let transmitted = if propagating {
        C64::from_re(2.0 * det::sqrt(incident) * det::sqrt(exit.re) * attenuation) / denominator
    } else {
        C64::ZERO
    };
    let power = reflected.norm_sq() + transmitted.norm_sq();
    if ![
        reflected.re,
        reflected.im,
        transmitted.re,
        transmitted.im,
        power,
    ]
    .iter()
    .all(|v| v.is_finite())
        || (power - 1.0).abs() > 256.0 * f64::EPSILON
    {
        return Err(DielectricError::InvalidInterface);
    }
    // Correct only bounded roundoff, jointly in the two Jones amplitudes.
    // This preserves relative phase and polarization positivity while closing
    // lossless work, including an interference null or a nearly opaque film.
    let scale = 1.0 / det::sqrt(power);
    Ok((reflected.scale(scale), transmitted.scale(scale)))
}

fn matrix(s: C64, p: C64) -> Mueller {
    let cross = s * p.conj();
    diagonal(s.norm_sq(), p.norm_sq(), cross.re, cross.im)
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

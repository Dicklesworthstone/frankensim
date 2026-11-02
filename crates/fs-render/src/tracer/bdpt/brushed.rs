//! Bidirectional transport for the native body-locked conductor finish.
//!
//! Reuse the camera owner's anisotropic GGX, visible-normal proposal, optical
//! tables and horizon support. Its geometric-measure radiance BSDF includes
//! the shading-normal projection and is not symmetric: forward light paths
//! need its adjoint. PBRT 3e, section 16.1.3, "Non-symmetry Due to Shading
//! Normals":
//! https://pbr-book.org/3ed-2018/Light_Transport_III_Bidirectional_Methods/The_Path-Space_Measurement_Equation

use super::*;

#[allow(clippy::too_many_arguments)]
pub(super) fn scattering(
    material: &Material,
    vertex: &Vertex<'_>,
    normal: Vec3,
    wo: Vec3,
    wi: Vec3,
    wavelength_nm: f64,
    mode: TransportMode,
) -> Result<f64, TracerError> {
    let incident_medium = vertex.medium_toward(wo).map(|entry| *entry.glass);
    // The radiance law is D G F / (4 ns.wo ng.wi); its importance
    // adjoint exchanges wo/wi, giving D G F / (4 ns.wi ng.wo).
    // Both directions are on the same geometric side for this opaque
    // reflection, so this exchange preserves the incident optical medium.
    let (outgoing, incoming) = match mode {
        TransportMode::Radiance => (wo, wi),
        TransportMode::Importance => (wi, wo),
    };
    opaque_bsdf_eval(
        material,
        normal,
        outgoing,
        incoming,
        wavelength_nm,
        incident_medium,
        vertex.brushed_frame.as_ref(),
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn sample(
    material: &Material,
    vertex: &Vertex<'_>,
    normal: Vec3,
    wo: Vec3,
    wavelength_nm: f64,
    mode: TransportMode,
    random: [f64; 4],
) -> Result<Option<ScatterSample>, TracerError> {
    let frame = vertex.brushed_frame.as_ref();
    let Some((direction, pdf_fwd)) = bsdf_sample(material, normal, wo, random[0], random[1], frame)
    else {
        return Ok(None);
    };
    let value = scattering(material, vertex, normal, wo, direction, wavelength_nm, mode)?;
    let weight = value * normal.dot(direction).max(0.0) / pdf_fwd;
    // The proposal is unchanged by the transport mode. Neither a shading-
    // normal adjoint nor rejected below-horizon samples renormalizes it.
    let pdf_rev = bsdf_pdf(material, normal, direction, wo, frame);
    Ok(
        (weight.is_finite() && weight > 0.0).then_some(ScatterSample {
            direction,
            weight,
            pdf_fwd_solid_angle: pdf_fwd,
            pdf_rev_solid_angle: pdf_rev,
            delta: false,
            event: None,
        }),
    )
}

#[cfg(test)]
mod tests;

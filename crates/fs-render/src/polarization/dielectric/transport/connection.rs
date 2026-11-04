//! Evaluated reciprocal surface connections for an intensity-measuring camera.
//!
//! The camera prefix retains its complete Stokes importance state. In opposite
//! travel frames sharing one horizontal axis, reciprocity contracts it with the
//! source through A = diag(1, 1, -1, 1), not a geometric frame reversal that
//! would also negate V. This is the same scalar reciprocal convention as the
//! native molecular connection, including coherent Fresnel phase.

use super::{BeamError, BeamMedium, PolarizationFrame, PolarizedBeam, unit};
use crate::dielectric::{
    DielectricCoatingStack, DielectricError, DielectricEvent, DielectricInterface, DielectricSurface,
};
use crate::polarization::{Mueller, Stokes};
use fs_geom::Vec3;
use fs_math::det;

/// A continuous radiance BSDF response and the original camera proposal.
///
/// `factor` is per steradian and includes neither input intensity, any source
/// cosine, source-leg transmission, nor camera index conversion. Multiply by
/// camera-prefix I, source irradiance on the actual geometric receiving plane,
/// and `(n_camera / n_camera_arrival)^2`. The rough dielectric's radiance law
/// already contains the local transmission eta factor: applying another
/// `(n_camera_arrival / n_light_arrival)^2` would count refraction twice.
///
/// `pdf` is the original scalar, polarization-independent camera solid-angle
/// proposal, including the native Fresnel branch probability and facet
/// rejection mass. It is not divided out of `factor`. Exact polarization
/// extinction can therefore have a positive PDF. Smooth/null delta transport
/// has zero continuous factor and PDF.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceConnectionResponse {
    pub factor: f64,
    pub pdf: f64,
    pub event: DielectricEvent,
}

impl SurfaceConnectionResponse {
    pub(crate) const fn zero(event: DielectricEvent) -> Self {
        Self { factor: 0.0, pdf: 0.0, event }
    }

    pub(crate) fn checked(factor: f64, pdf: f64, event: DielectricEvent) -> Result<Self, BeamError> {
        if !factor.is_finite() || factor < 0.0 || !pdf.is_finite() || pdf < 0.0 {
            return Err(BeamError::NumericalRange);
        }
        Ok(Self { factor, pdf, event })
    }
}

impl PolarizedBeam {
    /// Evaluate a rough dielectric at a fixed source-arrival direction.
    ///
    /// `self` is the original pre-scatter reciprocal camera beam. The source
    /// frame travels FROM the source TO this vertex; its Stokes intensity must
    /// be exactly one. `opposite_medium` is the OTHER original boundary medium
    /// even for reflection, where its index still determines Fresnel response.
    /// A scene consumer separately checks that the actual source arrives in
    /// the camera medium for reflection or the opposite medium for transmission.
    /// Retained gradient-index samples must be at this beam's physical origin.
    ///
    /// Reuses the native GGX evaluation, coherent stack ordering, Fresnel
    /// proposal and Mueller matrices. No direction or Fresnel event is sampled.
    /// The supplied normal is the original outward geometric normal, so exiting
    /// rays reverse authored coating order exactly as the native sampler does.
    pub fn evaluate_dielectric_connection(
        self,
        normal: Vec3,
        opposite_medium: impl Into<BeamMedium>,
        surface: DielectricSurface,
        source_frame: PolarizationFrame,
        normalized_source: Stokes,
    ) -> Result<SurfaceConnectionResponse, BeamError> {
        if surface.anisotropy().is_some() {
            return Err(DielectricError::UnsupportedAnisotropicTransport.into());
        }
        let reverse = normal.dot(self.frame.direction) > 0.0;
        self.evaluate_interface_connection(normal, opposite_medium.into(), surface.roughness_alpha(),
            surface.interface(self.wavelength_nm, reverse)?, source_frame, normalized_source)
    }

    /// Evaluate the same reciprocal rough response with an original borrowed
    /// mixed spectral stack. Its layers replace no hidden treatment: the
    /// supplied surface must be uncoated and declares only smooth/GGX roughness.
    /// Actual adjacent indices and authored outward-normal ordering are used
    /// exactly as in `sample_stack_interface`, including its original PDF.
    pub fn evaluate_stack_dielectric_connection(
        self,
        normal: Vec3,
        opposite_medium: impl Into<BeamMedium>,
        surface: DielectricSurface,
        stack: DielectricCoatingStack<'_>,
        source_frame: PolarizationFrame,
        normalized_source: Stokes,
    ) -> Result<SurfaceConnectionResponse, BeamError> {
        if surface.anisotropy().is_some() {
            return Err(DielectricError::UnsupportedAnisotropicTransport.into());
        }
        if surface.coating().is_some() {
            return Err(DielectricError::InvalidCoating.into());
        }
        let reverse = normal.dot(self.frame.direction) > 0.0;
        self.evaluate_interface_connection(normal, opposite_medium.into(), surface.roughness_alpha(),
            stack.interface(self.wavelength_nm, reverse)?, source_frame, normalized_source)
    }

    fn evaluate_interface_connection(
        self,
        normal: Vec3,
        opposite_medium: BeamMedium,
        roughness: Option<f64>,
        interface: DielectricInterface<'_>,
        source_frame: PolarizationFrame,
        normalized_source: Stokes,
    ) -> Result<SurfaceConnectionResponse, BeamError> {
        let camera = normalized_camera(self, normalized_source)?;
        let eta_i = camera.medium.phase_index_at(camera.origin, camera.wavelength_nm)?;
        let eta_t = opposite_medium.phase_index_at(camera.origin, camera.wavelength_nm)?;
        let normal = camera.incident_normal(normal)?;
        let wo = camera.frame.direction.scale(-1.0);
        let wi = source_frame.direction().scale(-1.0);
        let event = if normal.dot(wi) > 0.0 {
            DielectricEvent::Reflection
        } else {
            DielectricEvent::Transmission
        };
        let Some(alpha) = roughness else {
            return Ok(SurfaceConnectionResponse::zero(event));
        };
        let (evaluation, fresnel) = interface.evaluate_rough_polarized(normal, wo, wi, eta_i, eta_t, alpha)?;
        if evaluation.value == 0.0 {
            return SurfaceConnectionResponse::checked(0.0, evaluation.pdf, evaluation.event);
        }
        let ratio = if evaluation.event == DielectricEvent::Transmission { eta_t / eta_i } else { 1.0 };
        let mut facet = unit(Vec3::new(
            wo.x + ratio * wi.x, wo.y + ratio * wi.y, wo.z + ratio * wi.z,
        ))?;
        if facet.dot(normal) < 0.0 { facet = facet.scale(-1.0); }
        let fresnel = fresnel.ok_or(BeamError::NumericalRange)?;
        let matrix = match evaluation.event {
            DielectricEvent::Reflection => fresnel.reflection,
            DielectricEvent::Transmission => fresnel.transmission,
        };
        let mean = matrix.0[0][0];
        if !mean.is_finite() || mean <= 0.0 {
            return Err(BeamError::NumericalRange);
        }
        // Normalize the matrix first: dividing the BSDF by a very small
        // reflectance can overflow even when the final optical response fits.
        // The scalar evaluation already owns radiance eta-squared transport.
        let matrix = Mueller(matrix.0.map(|row| row.map(|value| value / mean)));
        let scattered = camera.scatter(facet, wi, matrix, evaluation.event, opposite_medium)?;
        let contraction = reciprocal_contraction(scattered, source_frame, normalized_source)?;
        let factor = checked_product(evaluation.value, contraction)?;
        SurfaceConnectionResponse::checked(factor, evaluation.pdf, evaluation.event)
    }
}

/// Admit a unit source and normalize the camera prefix before optical mixing.
/// The physical prefix intensity stays with the camera radiometric owner.
pub(crate) fn normalized_camera(
    camera: PolarizedBeam,
    normalized_source: Stokes,
) -> Result<PolarizedBeam, BeamError> {
    if normalized_source.i != 1.0 {
        return Err(BeamError::InvalidStokes);
    }
    normalized(normalized_source)?;
    let power = normalized(camera.power())?;
    PolarizedBeam::try_new(camera.ray().origin, camera.frame(), camera.wavelength_nm(), camera.medium(), power)
}

/// Contract an evaluated camera-side Mueller response with a unit arriving
/// source in the opposite travel frame. The camera response may be dark.
pub(crate) fn reciprocal_contraction(
    camera_response: PolarizedBeam,
    source_frame: PolarizationFrame,
    normalized_source: Stokes,
) -> Result<f64, BeamError> {
    if normalized_source.i != 1.0 {
        return Err(BeamError::InvalidStokes);
    }
    let input = normalized(normalized_source)?;
    let camera = camera_response.power();
    if camera.i == 0.0 { return Ok(0.0); }
    // The reconstructed facet can change the final direction at roundoff
    // scale. The caller's exact connection direction owns the output frame,
    // just as in the native rough conductor sampler.
    let reference = PolarizationFrame::try_new(
        source_frame.direction(), camera_response.frame().horizontal(),
    )?;
    let source = normalized(source_frame.express_in(input, reference))?;
    let response = normalized(camera)?;
    let a = [response.q, response.u, response.v];
    let b = [source.q, -source.u, source.v];
    // A positive-coherency overlap, factored to retain crossed-polarizer
    // extinction without cancelling two almost equal positive intensities.
    let sum_squares = a.into_iter().zip(b).map(|(a, b)| (a + b) * (a + b)).sum::<f64>();
    let overlap = 0.5 * (cone_slack(a)? + cone_slack(b)? + sum_squares);
    checked_product(camera.i, overlap)
}

pub(crate) fn checked_product(a: f64, b: f64) -> Result<f64, BeamError> {
    let value = a * b;
    if !a.is_finite() || a < 0.0 || !b.is_finite() || b < 0.0
        || !value.is_finite() || (a > 0.0 && b > 0.0 && value == 0.0) {
        return Err(BeamError::NumericalRange);
    }
    Ok(value)
}

fn normalized(input: Stokes) -> Result<Stokes, BeamError> {
    if input.i <= 0.0 || Stokes::try_new(input.i, input.q, input.u, input.v).is_none() {
        return Err(BeamError::InvalidStokes);
    }
    let mut result = Stokes { i: 1.0, q: input.q / input.i, u: input.u / input.i, v: input.v / input.i };
    let degree = det::hypot(det::hypot(result.q, result.u), result.v);
    if degree > 1.0 {
        result.q /= degree;
        result.u /= degree;
        result.v /= degree;
    }
    Ok(result)
}

fn cone_slack(polarization: [f64; 3]) -> Result<f64, BeamError> {
    let degree = det::hypot(det::hypot(polarization[0], polarization[1]), polarization[2]);
    if !degree.is_finite() || degree > 1.0 + 16.0 * f64::EPSILON {
        return Err(BeamError::NumericalRange);
    }
    Ok((1.0 - degree.min(1.0)) * (1.0 + degree.min(1.0)))
}

#[cfg(test)]
mod tests;

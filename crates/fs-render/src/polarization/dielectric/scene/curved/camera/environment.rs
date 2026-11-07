//! Explicit polarized environment radiance at the original escaped direction.
//!
//! RGB still owns only the declared scalar intensity reconstruction. A caller
//! supplies physical Stokes fractions, either shared by the complete map or
//! individually authored for each original texel. No polarization is inferred
//! from color, brightness, a sky name or the camera's polarization state.

use super::{
    CameraRadianceError, CurvedSceneError, EnvironmentRadiance, LAMBDA_MAX, LAMBDA_MIN,
    PathTermination, SourcePolarization, finite_nonnegative,
};
use crate::polarization::dielectric::transport::connection::{
    normalized_camera, reciprocal_contraction,
};
use crate::polarization::dielectric::transport::{PolarizationFrame, PolarizedBeam};
use fs_exec::Cx;
use fs_geom::Vec3;
use fs_math::det;

/// Directional source proposals with native visibility and complementary MIS.
pub mod sampling;

#[derive(Clone, Copy, Debug)]
enum Samples<'a> {
    Uniform(&'a SourcePolarization),
    Texels(&'a [SourcePolarization]),
}

/// The original map owns the direction-to-texel mapping. Trigonometry for its
/// rigid placement is prepared once and unused by ordinary scalar sources.
#[derive(Clone, Copy, Debug)]
pub(super) struct EnvironmentPolarization<'a> {
    samples: Samples<'a>,
    rotation: [f64; 2],
}

impl<'a> EnvironmentPolarization<'a> {
    fn new(samples: Samples<'a>, rotation_y_radians: f64) -> Self {
        let (sine, cosine) = det::sin_cos(rotation_y_radians);
        Self {
            samples,
            rotation: [sine, cosine],
        }
    }

    fn at(self, texel_index: usize) -> Option<&'a SourcePolarization> {
        match self.samples {
            Samples::Uniform(source) => Some(source),
            Samples::Texels(texels) => texels.get(texel_index),
        }
    }

    fn world_axis(self, local: Vec3) -> Vec3 {
        // EnvironmentMap uses phi_world = phi_local + rotation. Apply that
        // same proper rotation to the reference axis before ray projection.
        let [sine, cosine] = self.rotation;
        Vec3::new(
            cosine * local.x - sine * local.z,
            local.y,
            sine * local.x + cosine * local.z,
        )
    }

    fn response(self, texel_index: usize, beam: PolarizedBeam) -> Result<f64, CameraRadianceError> {
        let law = self
            .at(texel_index)
            .ok_or(CameraRadianceError::InvalidPolarization)?;
        let source = law.normalized_stokes(beam.wavelength_nm())?;
        if source.q == 0.0 && source.u == 0.0 && source.v == 0.0 {
            // An isotropic state has no reference-axis singularity. Preserve
            // the complete original scalar path, including every intensity bit.
            return Ok(1.0);
        }
        // Circular fractions have no linear orientation. A parallel declared
        // reference axis cannot make an otherwise physical circular sky dark
        // or undefined; use the already valid transverse transport axis.
        let axis = if source.q == 0.0 && source.u == 0.0 {
            beam.frame().horizontal()
        } else {
            self.world_axis(law.reference_axis())
        };
        let frame = PolarizationFrame::try_new(beam.ray().dir.scale(-1.0), axis)
            .map_err(CurvedSceneError::from)?;
        let camera = normalized_camera(beam, source).map_err(CurvedSceneError::from)?;
        reciprocal_contraction(camera, frame, source)
            .map_err(CurvedSceneError::from)
            .map_err(Into::into)
    }
}

fn admit_band(source: &SourcePolarization) -> Result<(), CameraRadianceError> {
    if source
        .wavelength_bounds_nm()
        .is_some_and(|(low, high)| low > LAMBDA_MIN || high < LAMBDA_MAX)
    {
        return Err(CameraRadianceError::InvalidPolarization);
    }
    Ok(())
}

impl<'a> EnvironmentRadiance<'a> {
    /// Attach an explicitly authored polarization law to every original texel.
    /// Its axis is in the unrotated map's Cartesian coordinates; its Stokes
    /// fractions are in the right-handed frame travelling FROM the environment
    /// toward the scene, opposite the escaped camera ray. In particular V uses
    /// that physical incoming-light convention, with no handedness guess.
    /// States with Q=U=0 have no linear orientation and need no projected axis;
    /// a singular declared axis refuses when a linear component is present.
    ///
    /// Spectral support must cover the complete visible RGB reconstruction.
    /// A query outside it refuses. This borrows the immutable law and allocates
    /// nothing; it introduces no atmospheric-scattering or calibration claim.
    pub fn with_polarization(
        mut self,
        source: &'a SourcePolarization,
    ) -> Result<Self, CameraRadianceError> {
        admit_band(source)?;
        self.polarization = Some(EnvironmentPolarization::new(
            Samples::Uniform(source),
            self.map.rotation_y_radians(),
        ));
        Ok(self)
    }

    /// Attach one independently authored spectral polarization law per
    /// original row-major texel. Length must equal width*height. Every law
    /// covers the full visible intensity band and uses map-local axes.
    ///
    /// Evaluation shares the exact original radiance lookup, including its
    /// seam, poles and map rotation. No neighboring polarization angles or
    /// Stokes coordinates in different frames are averaged. This is the same
    /// explicitly piecewise-constant angular model as the native RGB map;
    /// fractional Q/U/V retain each law's own convex spectral interpolation.
    /// Admission is cancellable, borrows storage and publishes only when all
    /// texels are valid. Camera queries allocate nothing and look up one law.
    pub fn with_polarization_map(
        mut self,
        cx: &Cx<'_>,
        texels: &'a [SourcePolarization],
    ) -> Result<Self, CameraRadianceError> {
        cx.checkpoint()?;
        let count = (self.map.width() as usize)
            .checked_mul(self.map.height() as usize)
            .ok_or(CameraRadianceError::InvalidPolarization)?;
        if texels.len() != count {
            return Err(CameraRadianceError::InvalidPolarization);
        }
        for (index, source) in texels.iter().enumerate() {
            if index % 32 == 0 {
                cx.checkpoint()?;
            }
            admit_band(source)?;
        }
        cx.checkpoint()?;
        self.polarization = Some(EnvironmentPolarization::new(
            Samples::Texels(texels),
            self.map.rotation_y_radians(),
        ));
        Ok(self)
    }

    /// The original authored law for one row-major texel, before map rotation.
    /// Uniform bindings return the same law for each valid texel; an absent
    /// declaration or an out-of-range texel returns None.
    pub fn polarization_at_texel(&self, texel_index: usize) -> Option<&'a SourcePolarization> {
        let count = (self.map.width() as usize).checked_mul(self.map.height() as usize)?;
        if texel_index >= count {
            return None;
        }
        self.polarization?.at(texel_index)
    }

    pub(super) fn evaluate_for_camera(
        &self,
        beam: PolarizedBeam,
    ) -> Result<(f64, f64), CameraRadianceError> {
        let Some(polarization) = self.polarization else {
            return Ok((self.evaluate(beam.ray().dir, beam.wavelength_nm())?, 1.0));
        };
        let (radiance, texel_index) = self.evaluate_texel(beam.ray().dir, beam.wavelength_nm())?;
        if radiance == 0.0 {
            return Ok((radiance, 1.0));
        }
        if beam.power().i == 0.0 {
            return Err(CameraRadianceError::Unresolved(
                PathTermination::Extinguished,
            ));
        }
        Ok((radiance, polarization.response(texel_index, beam)?))
    }
}

/// The geometric/native transport weight stays separately observable. Only
/// the endpoint source contraction changes its contribution to radiance.
pub(super) fn endpoint_radiance(
    emitted: f64,
    transport: f64,
    polarization: f64,
) -> Result<f64, CameraRadianceError> {
    if polarization == 0.0 {
        return Ok(0.0);
    }
    let direct = emitted * transport;
    if polarization == 1.0 {
        let result = finite_nonnegative(direct)?;
        if emitted > 0.0 && result == 0.0 {
            return Err(CameraRadianceError::NumericalRange);
        }
        return Ok(result);
    }
    if !polarization.is_finite() || polarization < 0.0 {
        return Err(CameraRadianceError::NumericalRange);
    }
    if emitted == 0.0 {
        return Ok(0.0);
    }
    let value = direct * polarization;
    if value.is_finite() && value > 0.0 {
        return Ok(value);
    }
    // A bounded analyzer response may rescue an overflowing intermediate.
    // All physical source/transport factors have already been admitted.
    let value = det::exp(det::ln(emitted) + det::ln(transport) + det::ln(polarization));
    if !value.is_finite() || value <= 0.0 {
        return Err(CameraRadianceError::NumericalRange);
    }
    Ok(value)
}

#[cfg(test)]
mod tests;

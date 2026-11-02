//! Add calibrated surface emission before the same native path scatters.
//!
//! The rendering equation adds Le at every encountered source, including a
//! surface that subsequently reflects or transmits. Its contribution uses the
//! arriving reciprocal Stokes I and that encounter's local index, never the
//! final path's index or the current surface's sampled BSDF/PDF.
//! See PBRT 4e, section 13.4:
//! <https://pbr-book.org/4ed/Light_Transport_I_Surface_Reflection/A_Better_Path_Tracer>.

use super::{
    CameraRadianceError, CurvedCameraScene, CurvedPathResult, CurvedSceneError, EmissionSides,
    PathTermination, finite_nonnegative, phase_index,
};

/// A stack-only accumulator. A callback failure stops the original walker and
/// retains the camera error without exposing camera policy in public power APIs.
pub(super) struct SurfaceEmission<'a, 's, 'm> {
    camera: &'a CurvedCameraScene<'s, 'm>,
    initial_index: f64,
    radiance: f64,
    error: Option<CameraRadianceError>,
}

impl<'a, 's, 'm> SurfaceEmission<'a, 's, 'm> {
    pub(super) fn new(camera: &'a CurvedCameraScene<'s, 'm>, initial_index: f64) -> Self {
        Self {
            camera,
            initial_index,
            radiance: 0.0,
            error: None,
        }
    }

    pub(super) fn observe(&mut self, path: &CurvedPathResult) -> Result<(), CurvedSceneError> {
        self.add_surface(path).map_err(|error| {
            self.error = Some(error);
            CurvedSceneError::NumericalRange
        })
    }

    fn add_surface(&mut self, path: &CurvedPathResult) -> Result<(), CameraRadianceError> {
        let contact = path
            .last_contact()
            .ok_or(CameraRadianceError::NumericalRange)?;
        let Some(emitter) = self
            .camera
            .emitters
            .get(contact.witness.boundary_id)
            .ok_or(CameraRadianceError::NumericalRange)?
        else {
            return Ok(());
        };
        if emitter.object_id != contact.object_id {
            return Err(CameraRadianceError::NumericalRange);
        }
        let beam = path.beam();
        // Sidedness belongs to the original geometry, not a face-forwarded
        // normal, sampled facet or authored normal perturbation.
        if emitter.sides == EmissionSides::Front
            && contact.witness.normal.dot(beam.ray().dir) >= 0.0
        {
            return Ok(());
        }
        let emitted = emitter.spectrum.evaluate(beam.wavelength_nm())?;
        if emitted == 0.0 {
            return Ok(());
        }
        if beam.power().i == 0.0 {
            return Err(CameraRadianceError::Unresolved(
                PathTermination::Extinguished,
            ));
        }
        let local_index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let ratio = self.initial_index / local_index;
        let weight = finite_nonnegative(beam.power().i * (ratio * ratio))?;
        let contribution = finite_nonnegative(emitted * weight)?;
        if contribution == 0.0 {
            return Err(CameraRadianceError::NumericalRange);
        }
        self.radiance = finite_nonnegative(self.radiance + contribution)?;
        Ok(())
    }

    pub(super) fn finish<T>(
        self,
        result: Result<T, CurvedSceneError>,
    ) -> Result<(T, f64), CameraRadianceError> {
        if let Some(error) = self.error {
            return Err(error);
        }
        Ok((result?, self.radiance))
    }
}

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
        self.observe_weighted(path, 1.0)
    }

    /// Apply a complementary camera-strategy weight to this emitter hit.
    /// Ordinary camera consumers retain weight one and their original arithmetic.
    pub(super) fn observe_weighted(
        &mut self,
        path: &CurvedPathResult,
        weight: f64,
    ) -> Result<(), CurvedSceneError> {
        self.add_surface(path, weight).map_err(|error| {
            self.error = Some(error);
            CurvedSceneError::NumericalRange
        })
    }

    fn add_surface(
        &mut self,
        path: &CurvedPathResult,
        mis_weight: f64,
    ) -> Result<(), CameraRadianceError> {
        if !mis_weight.is_finite() || !(0.0..=1.0).contains(&mis_weight) {
            return Err(CameraRadianceError::NumericalRange);
        }
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
        if emitted == 0.0 || mis_weight == 0.0 {
            return Ok(());
        }
        if beam.power().i == 0.0 {
            return Err(CameraRadianceError::Unresolved(
                PathTermination::Extinguished,
            ));
        }
        let local_index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let ratio = self.initial_index / local_index;
        let weight = beam.power().i * (ratio * ratio);
        let direct = (emitted * weight) * mis_weight;
        let contribution = if direct.is_finite() && direct > 0.0 {
            direct
        } else {
            // A small complementary weight can rescue an overflowing
            // intermediate, and the camera prefix can rescue underflow.
            let logarithm = fs_math::det::ln(emitted)
                + fs_math::det::ln(beam.power().i)
                + 2.0 * (fs_math::det::ln(self.initial_index) - fs_math::det::ln(local_index))
                + fs_math::det::ln(mis_weight);
            let value = fs_math::det::exp(logarithm);
            if !value.is_finite() || value <= 0.0 {
                return Err(CameraRadianceError::NumericalRange);
            }
            value
        };
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

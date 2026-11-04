//! Add surface and optional LTE volume emission before the native path scatters.
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

use super::super::CurvedTraceEvent;
use super::thermal::ThermalEmission;

#[cfg(feature = "volumes")]
use super::super::participating::{CurvedRayleighEvent, CurvedRayleighPath};

/// A stack-only accumulator. A callback failure stops the original walker and
/// retains the camera error without exposing camera policy in public power APIs.
pub(super) struct SurfaceEmission<'a, 's, 'm> {
    camera: &'a CurvedCameraScene<'s, 'm>,
    initial_index: f64,
    radiance: f64,
    error: Option<CameraRadianceError>,
    thermal: Option<ThermalEmission>,
    spatial_emission: f64,
}

impl<'a, 's, 'm> SurfaceEmission<'a, 's, 'm> {
    pub(super) fn new(camera: &'a CurvedCameraScene<'s, 'm>, initial_index: f64) -> Self {
        Self {
            camera,
            initial_index,
            radiance: 0.0,
            error: None,
            thermal: (!camera.volume_emitters.is_empty()).then(ThermalEmission::new),
            spatial_emission: 0.0,
        }
    }

    pub(super) fn observe(&mut self, path: &CurvedPathResult) -> Result<(), CurvedSceneError> {
        self.observe_weighted(path, 1.0)
    }

    pub(super) fn observe_event(
        &mut self,
        event: CurvedTraceEvent<'_>,
    ) -> Result<(), CurvedSceneError> {
        match event {
            CurvedTraceEvent::Arrival(path) => self.observe(path),
            CurvedTraceEvent::Departure { path, .. } => {
                self.depart(path);
                Ok(())
            }
        }
    }

    #[cfg(feature = "volumes")]
    pub(super) fn observe_rayleigh_event(
        &mut self,
        event: CurvedRayleighEvent<'_>,
    ) -> Result<(), CurvedSceneError> {
        match event {
            CurvedRayleighEvent::Surface(event) => self.observe_event(event),
            CurvedRayleighEvent::VolumeArrival { path, .. } => self.observe_volume(path),
            CurvedRayleighEvent::VolumeDeparture { path, .. } => {
                self.depart(path);
                Ok(())
            }
        }
    }

    /// The native response has already applied its Mueller/PDF weight.
    pub(super) fn depart(&mut self, path: &CurvedPathResult) {
        if let Some(thermal) = &mut self.thermal {
            thermal.depart(path);
        }
    }

    #[cfg(feature = "volumes")]
    pub(super) fn observe_volume(
        &mut self,
        path: &CurvedPathResult,
    ) -> Result<(), CurvedSceneError> {
        self.add_thermal(path, false).map_err(|error| {
            self.error = Some(error);
            CurvedSceneError::NumericalRange
        })
    }

    fn add_thermal(
        &mut self,
        path: &CurvedPathResult,
        terminal: bool,
    ) -> Result<(), CameraRadianceError> {
        if !self.camera.spatial_volume_emitters.is_empty() {
            let increment =
                finite_nonnegative(path.reduced_volume_emission - self.spatial_emission)?;
            if increment != 0.0 {
                let emitted =
                    super::thermal::weighted_source(increment, self.initial_index, 1.0, 1.0)?;
                self.radiance = finite_nonnegative(self.radiance + emitted)?;
            }
            self.spatial_emission = path.reduced_volume_emission;
        }
        if let Some(thermal) = &mut self.thermal {
            let emitted = thermal.arrive(self.camera, self.initial_index, path, terminal)?;
            // Avoid changing even a retained zero's sign when no source emits.
            if emitted != 0.0 {
                self.radiance = finite_nonnegative(self.radiance + emitted)?;
            }
        }
        Ok(())
    }

    /// Apply a complementary camera-strategy weight to this emitter hit.
    /// Ordinary camera consumers retain weight one and their original arithmetic.
    pub(super) fn observe_weighted(
        &mut self,
        path: &CurvedPathResult,
        weight: f64,
    ) -> Result<(), CurvedSceneError> {
        // Area-source MIS applies only to the surface term. Thermal emission
        // has no competing source proposal and retains its full path weight.
        self.add_thermal(path, false)
            .and_then(|()| self.add_surface(path, weight))
            .map_err(|error| {
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
        let emitter = self
            .camera
            .emitters
            .get(contact.witness.boundary_id)
            .ok_or(CameraRadianceError::NumericalRange)?;
        let thermal = self
            .camera
            .diffuse_thermal_emitters
            .get(contact.witness.boundary_id)
            .and_then(Option::as_ref);
        let (object_id, sides) = match (emitter, thermal) {
            (Some(emitter), None) => (emitter.object_id, emitter.sides),
            (None, Some(thermal)) => (thermal.object_id(), EmissionSides::Both),
            (None, None) => return Ok(()),
            (Some(_), Some(_)) => return Err(CameraRadianceError::NumericalRange),
        };
        if object_id != contact.object_id {
            return Err(CameraRadianceError::NumericalRange);
        }
        let beam = path.beam();
        // Sidedness belongs to the original geometry, not a face-forwarded
        // normal, sampled facet or authored normal perturbation.
        if sides == EmissionSides::Front && contact.witness.normal.dot(beam.ray().dir) >= 0.0 {
            return Ok(());
        }
        let emitted = match (emitter, thermal) {
            (Some(emitter), None) => emitter.spectrum.evaluate(beam.wavelength_nm())?,
            (None, Some(thermal)) => {
                // This source's original Kirchhoff/Planck law emits into
                // vacuum. Keep actual local phase-index admission separate
                // from the reciprocal camera-to-source radiance conversion.
                if phase_index(beam.medium(), beam.wavelength_nm())? != 1.0 {
                    return Err(CameraRadianceError::UnsupportedSurface { object_id });
                }
                thermal
                    .surface()
                    .radiance_w_m2_sr_nm(beam.wavelength_nm())
                    .map_err(CameraRadianceError::Thermal)?
            }
            _ => return Err(CameraRadianceError::NumericalRange),
        };
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

    pub(super) fn finish(
        self,
        result: Result<CurvedPathResult, CurvedSceneError>,
    ) -> Result<(CurvedPathResult, f64), CameraRadianceError> {
        self.finish_with(result, |path| path)
    }

    #[cfg(feature = "volumes")]
    pub(super) fn finish_scattering(
        self,
        result: Result<CurvedRayleighPath, CurvedSceneError>,
    ) -> Result<(CurvedRayleighPath, f64), CameraRadianceError> {
        self.finish_with(result, CurvedRayleighPath::path)
    }

    fn finish_with<T>(
        mut self,
        result: Result<T, CurvedSceneError>,
        path: fn(&T) -> &CurvedPathResult,
    ) -> Result<(T, f64), CameraRadianceError> {
        if let Some(error) = self.error.take() {
            return Err(error);
        }
        let result = result?;
        let terminal = path(&result);
        if matches!(
            terminal.termination(),
            PathTermination::Opaque
                | PathTermination::Escaped
                | PathTermination::RejectedMicrofacet
        ) {
            self.add_thermal(terminal, true)?;
        }
        Ok((result, self.radiance))
    }
}

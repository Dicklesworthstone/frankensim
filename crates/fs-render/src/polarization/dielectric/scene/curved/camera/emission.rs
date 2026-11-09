//! Add surface and optional LTE volume emission before the native path scatters.
//!
//! The rendering equation adds Le at every encountered source, including a
//! surface that subsequently reflects or transmits. Its contribution uses the
//! arriving reciprocal Stokes response and that encounter's local index, never the
//! final path's index or the current surface's sampled BSDF/PDF.
//! See PBRT 4e, section 13.4:
//! <https://pbr-book.org/4ed/Light_Transport_I_Surface_Reflection/A_Better_Path_Tracer>.

use super::{
    CameraRadianceError, CurvedCameraScene, CurvedPathResult, CurvedSceneError, EmissionSides,
    PathTermination, finite_nonnegative, phase_index,
};

use super::super::CurvedTraceEvent;
use super::thermal::ThermalEmission;
use crate::polarization::dielectric::transport::{PolarizationFrame, PolarizedBeam};
use crate::polarization::dielectric::transport::connection::{
    normalized_camera, reciprocal_contraction,
};

#[cfg(feature = "volumes")]
use super::super::participating::{CurvedRayleighEvent, CurvedRayleighPath};

/// Positive incoherent s/p sources must remain separate through a crossed
/// analyzer: forming I=s+p and Q=s-p first can erase the weaker channel.
/// The result is relative to total emitted radiance; physical camera prefix I
/// remains in the original radiometric/logarithmic accumulation below.
pub(super) fn thermal_polarization_weight(
    beam: PolarizedBeam,
    frame: PolarizationFrame,
    radiance: [f64; 2],
) -> Result<f64, CameraRadianceError> {
    use crate::polarization::Stokes;
    let total = radiance[0] + radiance[1];
    if !total.is_finite()
        || total <= 0.0
        || radiance.iter().any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(CameraRadianceError::NumericalRange);
    }
    if radiance[0] == radiance[1] {
        return Ok(1.0);
    }
    let camera = normalized_camera(beam, Stokes::unpolarized(1.0))
        .map_err(CurvedSceneError::from)?;
    let s = reciprocal_contraction(camera, frame,
        Stokes { i: 1.0, q: 1.0, u: 0.0, v: 0.0 })
        .map_err(CurvedSceneError::from)?;
    let p = reciprocal_contraction(camera, frame,
        Stokes { i: 1.0, q: -1.0, u: 0.0, v: 0.0 })
        .map_err(CurvedSceneError::from)?;
    let weight = (radiance[0] / total) * s + (radiance[1] / total) * p;
    if !weight.is_finite()
        || weight < 0.0
        || (weight == 0.0
            && ((radiance[0] > 0.0 && s > 0.0) || (radiance[1] > 0.0 && p > 0.0)))
    {
        return Err(CameraRadianceError::NumericalRange);
    }
    Ok(weight)
}

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

    /// An explicit source owner supplies already weighted physical radiance.
    /// Preserve its failure until `finish`; never return a partial estimate.
    pub(super) fn add_radiance(
        &mut self,
        contribution: Result<f64, CameraRadianceError>,
    ) -> Result<(), CurvedSceneError> {
        contribution.and_then(finite_nonnegative).and_then(|value| {
            if value != 0.0 {
                self.radiance = finite_nonnegative(self.radiance + value)?;
            }
            Ok(())
        }).map_err(|error| {
            self.error = Some(error);
            CurvedSceneError::NumericalRange
        })
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
        // Both native spatial-field LTE and finite absorption-cell LTE feed
        // this reduced-radiance counter, already weighted by the actual full
        // reciprocal Mueller prefix. Consume each increment exactly once.
        if path.reduced_volume_emission != self.spatial_emission {
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
        self.observe_weighted_modulated(path, weight, 1.0)
    }

    /// A finite source's prescribed radiance history changes only its own Le.
    /// Steady LTE emission and the original reciprocal continuation are unchanged.
    pub(super) fn observe_modulated(
        &mut self,
        path: &CurvedPathResult,
        source_scale: f64,
    ) -> Result<(), CurvedSceneError> {
        self.observe_weighted_modulated(path, 1.0, source_scale)
    }

    /// Combine the original emitter-hit MIS weight with that source's own
    /// retarded history before range recovery. Neither factor changes LTE.
    pub(super) fn observe_weighted_modulated(
        &mut self,
        path: &CurvedPathResult,
        weight: f64,
        source_scale: f64,
    ) -> Result<(), CurvedSceneError> {
        self.add_thermal(path, false)
            .and_then(|()| self.add_surface(path, weight, source_scale))
            .map_err(|error| {
                self.error = Some(error);
                CurvedSceneError::NumericalRange
            })
    }

    fn add_surface(
        &mut self,
        path: &CurvedPathResult,
        mis_weight: f64,
        source_scale: f64,
    ) -> Result<(), CameraRadianceError> {
        if !mis_weight.is_finite()
            || !(0.0..=1.0).contains(&mis_weight)
            || !source_scale.is_finite()
            || source_scale < 0.0
        {
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
        let beam = path.beam();
        let coating = self.camera.coating_thermal_emitters
            .get(contact.witness.boundary_id)
            .and_then(|faces| {
                let side = usize::from(contact.witness.normal.dot(beam.ray().dir) > 0.0);
                faces[side].as_ref()
            });
        let (object_id, sides) = match (emitter, thermal, coating) {
            (Some(emitter), None, None) => (emitter.object_id, emitter.sides),
            (None, Some(thermal), None) => (thermal.object_id(), EmissionSides::Both),
            (None, None, Some(thermal)) => (thermal.object_id(), EmissionSides::Both),
            (None, None, None) => return Ok(()),
            _ => return Err(CameraRadianceError::NumericalRange),
        };
        if object_id != contact.object_id {
            return Err(CameraRadianceError::NumericalRange);
        }
        // Sidedness belongs to the original geometry, not a face-forwarded
        // normal, sampled facet or authored normal perturbation.
        if sides == EmissionSides::Front && contact.witness.normal.dot(beam.ray().dir) >= 0.0 {
            return Ok(());
        }
        let mut directional_state = None;
        let emitted = match (emitter, thermal, coating) {
            (Some(emitter), None, None) => emitter.spectrum.evaluate(beam.wavelength_nm())?,
            (None, Some(thermal), None) => {
                // This source's original Kirchhoff/Planck law emits into
                // vacuum. Keep actual local phase-index admission separate
                // from the reciprocal camera-to-source radiance conversion.
                if phase_index(beam.medium(), beam.wavelength_nm())? != 1.0 {
                    return Err(CameraRadianceError::UnsupportedSurface { object_id });
                }
                thermal
                    .radiance_along(
                        beam.wavelength_nm(),
                        contact.witness.normal,
                        beam.ray().dir.scale(-1.0),
                    )
                    .map_err(CameraRadianceError::Thermal)?
            }
            (None, None, Some(thermal)) => {
                if phase_index(beam.medium(), beam.wavelength_nm())? != 1.0 {
                    return Err(CameraRadianceError::UnsupportedSurface { object_id });
                }
                let (radiance, frame) = thermal.sample(
                    beam.wavelength_nm(), beam.ray().dir.scale(-1.0),
                    contact.witness.normal, beam.frame().horizontal(),
                )?;
                directional_state = Some((radiance, frame));
                radiance[0] + radiance[1]
            }
            _ => return Err(CameraRadianceError::NumericalRange),
        };
        if emitted == 0.0 || mis_weight == 0.0 || source_scale == 0.0 {
            return Ok(());
        }
        if beam.power().i == 0.0 {
            return Err(CameraRadianceError::Unresolved(
                PathTermination::Extinguished,
            ));
        }
        let source_polarization = if let Some(polarization) =
            emitter.and_then(|emitter| emitter.polarization())
        {
            let source = polarization.normalized_stokes(beam.wavelength_nm())?;
            if source.q == 0.0 && source.u == 0.0 && source.v == 0.0 {
                // Retain every legacy intensity bit and its arbitrary frame.
                None
            } else {
                let axis = self.camera.emitter_polarization_axes
                    .get(contact.witness.boundary_id)
                    .and_then(|axis| *axis)
                    .ok_or(CameraRadianceError::NumericalRange)?;
                let frame = PolarizationFrame::try_new(beam.ray().dir.scale(-1.0), axis)
                    .map_err(CurvedSceneError::from)?;
                Some((source, frame))
            }
        } else {
            None
        };
        let polarization_weight = if let Some((radiance, frame)) = directional_state {
            thermal_polarization_weight(*beam, frame, radiance)?
        } else if let Some((source, frame)) = source_polarization {
            // Keep physical prefix I outside the bounded polarization
            // contraction so the existing logarithmic recovery can rescue
            // extreme source/prefix products. Reciprocity negates U only.
            let camera = normalized_camera(*beam, source).map_err(CurvedSceneError::from)?;
            reciprocal_contraction(camera, frame, source).map_err(CurvedSceneError::from)?
        } else {
            1.0
        };
        // A crossed source/analyzer state is exact zero emitted contribution,
        // not an extinguished path. The original material still scatters.
        if polarization_weight == 0.0 {
            return Ok(());
        }
        let local_index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let ratio = self.initial_index / local_index;
        let weight = beam.power().i * (ratio * ratio);
        let direct = (emitted * weight) * mis_weight;
        let direct = if polarization_weight == 1.0 {
            direct
        } else {
            direct * polarization_weight
        };
        // Keep the complete static arithmetic, including signed zero and the
        // original overflow recovery, when the history is exactly steady.
        let direct = if source_scale == 1.0 {
            direct
        } else {
            direct * source_scale
        };
        let contribution = if direct.is_finite() && direct > 0.0 {
            direct
        } else {
            // A small complementary weight can rescue an overflowing
            // intermediate, and the camera prefix can rescue underflow.
            let logarithm = fs_math::det::ln(emitted)
                + fs_math::det::ln(beam.power().i)
                + 2.0 * (fs_math::det::ln(self.initial_index) - fs_math::det::ln(local_index))
                + fs_math::det::ln(mis_weight);
            let logarithm = if polarization_weight == 1.0 {
                logarithm
            } else {
                logarithm + fs_math::det::ln(polarization_weight)
            };
            let logarithm = if source_scale == 1.0 {
                logarithm
            } else {
                logarithm + fs_math::det::ln(source_scale)
            };
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

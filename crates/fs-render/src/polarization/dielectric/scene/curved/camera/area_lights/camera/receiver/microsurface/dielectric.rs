//! Full two-sided glass connections at the original native receiving contact.
//! The source and camera may arrive in different adjacent media. The core
//! response owns the net local radiance eta factor; the camera-prefix conversion
//! remains with the original radiometric accumulator.

use super::*;
use crate::dielectric::DielectricSurface;
use crate::polarization::dielectric::transport::BeamMedium;
use crate::polarization::dielectric::transport::microsurface::dielectric::SmithDielectricError;

impl SurfaceReceiver {
    /// Identify the full-glass continuation law even when its continuous
    /// connection estimator is unavailable (equal-index mixed measures).
    pub(in super::super::super::super::super) fn uses_smith_dielectric(
        self,
        camera: &CurvedCameraScene<'_, '_>,
    ) -> Result<bool, AreaLightCameraError> {
        Ok(self.smith_dielectric_media(camera)?.is_some())
    }

    /// Original physical adjacent IDs, with no random draws or optical solve.
    pub(in super::super::super::super::super) fn smith_dielectric_media(
        self,
        camera: &CurvedCameraScene<'_, '_>,
    ) -> Result<Option<[usize; 2]>, AreaLightCameraError> {
        let native = camera
            .scene
            .bindings
            .get(self.contact.witness.boundary_id)
            .ok_or(AreaLightCameraError::NumericalRange)?;
        if native.object_id != self.contact.object_id {
            return Err(AreaLightCameraError::NumericalRange);
        }
        Ok(match native.interaction {
            Interaction::Dielectric {
                outside,
                inside,
                surface,
            } if !surface.is_delta() => Some([outside, inside]),
            _ => None,
        })
    }

    pub(in super::super::super::super::super) fn directional_source_medium(
        self,
        camera: &CurvedCameraScene<'_, '_>,
        direction: Vec3,
    ) -> Result<usize, AreaLightCameraError> {
        let native = camera
            .scene
            .bindings
            .get(self.contact.witness.boundary_id)
            .ok_or(AreaLightCameraError::NumericalRange)?;
        if native.object_id != self.contact.object_id {
            return Err(AreaLightCameraError::NumericalRange);
        }
        arrival_medium(native.interaction, self.contact, direction)
    }

    /// Equal-index rough glass retains its complete path-only law. A coated
    /// equal-index interface has a straight-through atom as well as continuous
    /// reflection, so callers must not replace its endpoint with continuous NEE.
    pub(in super::super::super::super::super) fn supports_smith_dielectric_connections(
        self,
        camera: &CurvedCameraScene<'_, '_>,
    ) -> Result<bool, AreaLightCameraError> {
        let Some((opposite, surface)) =
            self.smith_dielectric_state(camera, self.incoming.ray().origin)?
        else {
            return Ok(false);
        };
        self.incoming
            .smith_dielectric_connection_is_continuous(opposite, surface)
            .map_err(CurvedSceneError::from)
            .map_err(CameraRadianceError::from)
            .map_err(Into::into)
    }

    #[allow(clippy::too_many_arguments)]
    pub(in super::super::super::super::super) fn smith_dielectric_response(
        self,
        cx: &Cx<'_>,
        camera: &CurvedCameraScene<'_, '_>,
        connection: &BundleConnection,
        binding: usize,
        settings: SmithWalkSettings,
        samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<(Vec3, f64), AreaLightCameraError> {
        cx.checkpoint()?;
        let path = connection.path();
        let contact = *path
            .last_contact()
            .ok_or(AreaLightCameraError::PathMismatch { binding })?;
        if !path.is_surface_receiver() || !same_contact(self.contact, contact) {
            return Err(AreaLightCameraError::PathMismatch { binding });
        }
        let bundle = path.bundle();
        let endpoint = bundle.propagation().propagation();
        if self.directional_source_medium(camera, endpoint.ray.dir)? != path.medium_index() {
            return Err(AreaLightCameraError::PathMismatch { binding });
        }
        let frame = PolarizationFrame::try_new(endpoint.ray.dir, endpoint.transverse_axis)
            .map_err(CurvedSceneError::from)
            .map_err(CameraRadianceError::from)?;
        let factor = self.smith_dielectric_response_at(
            cx,
            camera,
            endpoint.ray.origin,
            frame,
            bundle.interface_state(),
            settings,
            samples,
        )?;
        Ok((face_normal(contact, endpoint.ray.dir), factor))
    }

    #[allow(clippy::too_many_arguments)]
    pub(in super::super::super::super::super) fn directional_smith_dielectric_response(
        self,
        cx: &Cx<'_>,
        camera: &CurvedCameraScene<'_, '_>,
        source_frame: PolarizationFrame,
        source_stokes: Stokes,
        settings: SmithWalkSettings,
        samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<(Vec3, f64, usize), AreaLightCameraError> {
        let medium = self.directional_source_medium(camera, source_frame.direction())?;
        let factor = self.smith_dielectric_response_at(
            cx,
            camera,
            self.incoming.ray().origin,
            source_frame,
            source_stokes,
            settings,
            samples,
        )?;
        Ok((
            face_normal(self.contact, source_frame.direction()),
            factor,
            medium,
        ))
    }

    /// Bind physical anisotropic body axes once, retaining the actual phase
    /// indices and original camera-side medium at this geometric point.
    fn smith_dielectric_state(
        self,
        camera: &CurvedCameraScene<'_, '_>,
        position: Point3,
    ) -> Result<Option<(BeamMedium, DielectricSurface)>, AreaLightCameraError> {
        let native = camera
            .scene
            .bindings
            .get(self.contact.witness.boundary_id)
            .ok_or(AreaLightCameraError::NumericalRange)?;
        if native.object_id != self.contact.object_id
            || arrival_medium(native.interaction, self.contact, self.incoming.ray().dir)?
                != self.medium
        {
            return Err(AreaLightCameraError::NumericalRange);
        }
        let Interaction::Dielectric {
            outside,
            inside,
            surface,
        } = native.interaction
        else {
            return Ok(None);
        };
        if surface.is_delta() {
            return Ok(None);
        }
        let opposite = if self.medium == outside {
            inside
        } else {
            outside
        };
        let medium = camera
            .scene
            .medium_at(opposite, position, self.incoming.wavelength_nm())
            .map_err(CameraRadianceError::from)?;
        let instance = camera
            .scene
            .geometry
            .instances()
            .get(self.contact.witness.boundary_id)
            .ok_or(AreaLightCameraError::NumericalRange)?;
        let surface = bind_dielectric_surface(surface, instance.transform())
            .map_err(CurvedSceneError::from)
            .map_err(CameraRadianceError::from)?;
        Ok(Some((medium, surface)))
    }

    #[allow(clippy::too_many_arguments)]
    fn smith_dielectric_response_at(
        self,
        cx: &Cx<'_>,
        camera: &CurvedCameraScene<'_, '_>,
        position: Point3,
        source_frame: PolarizationFrame,
        source_stokes: Stokes,
        settings: SmithWalkSettings,
        samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<f64, AreaLightCameraError> {
        cx.checkpoint()?;
        let (opposite, surface) = self
            .smith_dielectric_state(camera, position)?
            .ok_or(AreaLightCameraError::NumericalRange)?;
        let incoming = camera
            .scene
            .launch_beam(
                position,
                self.incoming.frame(),
                self.incoming.wavelength_nm(),
                self.medium,
                self.incoming.power(),
            )
            .map_err(CameraRadianceError::from)?;
        let response = if let Some(stack) = camera
            .scene
            .coating_stacks
            .get(self.contact.witness.boundary_id)
            .copied()
            .flatten()
        {
            incoming.evaluate_smith_stack_dielectric_connection(
                cx,
                self.contact.witness.normal,
                opposite,
                surface,
                stack,
                source_frame,
                source_stokes,
                settings,
                samples,
            )
        } else {
            incoming.evaluate_smith_dielectric_connection(
                cx,
                self.contact.witness.normal,
                opposite,
                surface,
                source_frame,
                source_stokes,
                settings,
                samples,
            )
        };
        response
            .map(|response| response.factor_per_sr())
            .map_err(smith_dielectric_error)
    }
}

fn smith_dielectric_error(error: SmithDielectricError) -> AreaLightCameraError {
    CameraRadianceError::from(CurvedSceneError::from(SceneError::from(error))).into()
}

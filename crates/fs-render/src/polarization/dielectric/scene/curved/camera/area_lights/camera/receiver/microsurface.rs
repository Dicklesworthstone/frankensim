//! Shared full Smith receivers for finite and infinite-direction sources.
//! No complete-walk directional PDF is exposed.

use super::*;
use crate::polarization::dielectric::transport::microsurface::connection::SmithConnectionResponse;
use crate::polarization::dielectric::transport::microsurface::{SmithBeamError, SmithWalkSettings};
use fs_exec::Cx;

mod dielectric;

impl SurfaceReceiver {
    #[allow(clippy::too_many_arguments)]
    pub(in super::super::super::super) fn smith_response(
        self,
        cx: &Cx<'_>,
        camera: &CurvedCameraScene<'_, '_>,
        connection: &BundleConnection,
        binding: usize,
        settings: SmithWalkSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<(Vec3, f64), AreaLightCameraError> {
        let native = camera
            .scene
            .bindings
            .get(self.contact.witness.boundary_id)
            .ok_or(AreaLightCameraError::PathMismatch { binding })?;
        if !matches!(
            native.interaction,
            Interaction::RoughConductor { .. } | Interaction::BrushedConductor { .. }
        ) {
            let response = self.response(camera, connection, binding)?;
            return Ok((response.normal, response.factor));
        }
        cx.checkpoint()?;
        let path = connection.path();
        let contact = *path
            .last_contact()
            .ok_or(AreaLightCameraError::PathMismatch { binding })?;
        let bundle = path.bundle();
        let endpoint = bundle.propagation().propagation();
        if !path.is_surface_receiver()
            || !same_contact(self.contact, contact)
            || native.object_id != contact.object_id
            || arrival_medium(native.interaction, contact, self.incoming.ray().dir)? != self.medium
            || arrival_medium(native.interaction, contact, endpoint.ray.dir)? != path.medium_index()
        {
            return Err(AreaLightCameraError::PathMismatch { binding });
        }
        let source_frame = PolarizationFrame::try_new(endpoint.ray.dir, endpoint.transverse_axis)
            .map_err(CurvedSceneError::from)
            .map_err(CameraRadianceError::from)?;
        let response = self.smith_response_at(
            cx,
            camera,
            endpoint.ray.origin,
            source_frame,
            bundle.interface_state(),
            settings,
            samples,
        )?;
        Ok((
            face_normal(contact, endpoint.ray.dir),
            response.factor_per_sr(),
        ))
    }
    pub(in super::super::super::super) fn uses_smith(
        self,
        camera: &CurvedCameraScene<'_, '_>,
    ) -> Result<bool, AreaLightCameraError> {
        let native = camera
            .scene
            .bindings
            .get(self.contact.witness.boundary_id)
            .ok_or(AreaLightCameraError::NumericalRange)?;
        if native.object_id != self.contact.object_id {
            return Err(AreaLightCameraError::NumericalRange);
        }
        Ok(match native.interaction {
            Interaction::RoughConductor { surface, .. } => !surface.is_delta(),
            Interaction::BrushedConductor { .. } => true,
            _ => false,
        })
    }

    /// Infinite-direction source response and its exact adjacent medium.
    /// The returned first-order geometric density is distinct from the full,
    /// stochastic BRDF and is never presented as its outgoing proposal PDF.
    #[allow(clippy::too_many_arguments)]
    pub(in super::super::super::super) fn directional_smith_response(
        self,
        cx: &Cx<'_>,
        camera: &CurvedCameraScene<'_, '_>,
        source_frame: PolarizationFrame,
        source_stokes: Stokes,
        settings: SmithWalkSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<(Vec3, SmithConnectionResponse, usize), AreaLightCameraError> {
        let native = camera
            .scene
            .bindings
            .get(self.contact.witness.boundary_id)
            .ok_or(AreaLightCameraError::NumericalRange)?;
        let medium = arrival_medium(native.interaction, self.contact, source_frame.direction())?;
        let response = self.smith_response_at(
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
            response,
            medium,
        ))
    }

    /// Geometry-only endpoint partition query. No random walk or Fresnel solve.
    pub(in super::super::super::super) fn smith_first_reflection_density(
        self,
        camera: &CurvedCameraScene<'_, '_>,
        outgoing: Vec3,
    ) -> Result<f64, AreaLightCameraError> {
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
        match native.interaction {
            Interaction::RoughConductor { surface, .. } => {
                self.incoming.smith_conductor_first_reflection_density(
                    self.contact.witness.normal,
                    surface,
                    outgoing,
                )
            }
            Interaction::BrushedConductor { detail, .. } => {
                let (transform, local) = self.smith_placement(camera)?;
                self.incoming
                    .instanced_smith_brushed_first_reflection_density(
                        &detail,
                        transform,
                        local,
                        self.contact.witness.normal,
                        outgoing,
                    )
            }
            _ => return Err(AreaLightCameraError::NumericalRange),
        }
        .map_err(smith_error)
    }

    fn smith_placement(
        self,
        camera: &CurvedCameraScene<'_, '_>,
    ) -> Result<(crate::instances::RigidTransform, Point3), AreaLightCameraError> {
        let instance = camera
            .scene
            .geometry
            .instances()
            .get(self.contact.witness.boundary_id)
            .ok_or(AreaLightCameraError::NumericalRange)?;
        let transform = instance.transform();
        let local = transform
            .inverse()
            .map_err(SceneError::from)
            .map_err(CurvedSceneError::from)
            .map_err(CameraRadianceError::from)?
            .transform_point(self.contact.witness.surface_point);
        Ok((transform, local))
    }

    #[allow(clippy::too_many_arguments)]
    fn smith_response_at(
        self,
        cx: &Cx<'_>,
        camera: &CurvedCameraScene<'_, '_>,
        position: Point3,
        source_frame: PolarizationFrame,
        source_stokes: Stokes,
        settings: SmithWalkSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<SmithConnectionResponse, AreaLightCameraError> {
        cx.checkpoint()?;
        let native = camera
            .scene
            .bindings
            .get(self.contact.witness.boundary_id)
            .ok_or(AreaLightCameraError::NumericalRange)?;
        if native.object_id != self.contact.object_id
            || arrival_medium(native.interaction, self.contact, self.incoming.ray().dir)?
                != self.medium
            || arrival_medium(native.interaction, self.contact, source_frame.direction())?
                != self.medium
        {
            return Err(AreaLightCameraError::NumericalRange);
        }
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
        match native.interaction {
            Interaction::RoughConductor {
                optics, surface, ..
            } => incoming.evaluate_smith_conductor_connection(
                cx,
                self.contact.witness.normal,
                &camera.scene.conductors[optics],
                surface,
                source_frame,
                source_stokes,
                settings,
                samples,
            ),
            Interaction::BrushedConductor { optics, detail, .. } => {
                let (transform, local) = self.smith_placement(camera)?;
                incoming.evaluate_instanced_smith_brushed_conductor_connection(
                    cx,
                    &camera.scene.conductors[optics],
                    &detail,
                    transform,
                    local,
                    self.contact.witness.normal,
                    source_frame,
                    source_stokes,
                    settings,
                    samples,
                )
            }
            _ => return Err(AreaLightCameraError::NumericalRange),
        }
        .map_err(smith_error)
    }
}

fn smith_error(error: SmithBeamError) -> AreaLightCameraError {
    CameraRadianceError::from(CurvedSceneError::from(SceneError::from(error))).into()
}

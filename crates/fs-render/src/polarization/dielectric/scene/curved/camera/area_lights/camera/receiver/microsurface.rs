//! Complete Smith responses for point fixtures; no directional PDF is exposed.

use super::*;
use crate::polarization::dielectric::transport::microsurface::{SmithBeamError, SmithWalkSettings};
use fs_exec::Cx;

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
        let incoming = camera
            .scene
            .launch_beam(
                endpoint.ray.origin,
                self.incoming.frame(),
                self.incoming.wavelength_nm(),
                self.medium,
                self.incoming.power(),
            )
            .map_err(CameraRadianceError::from)?;
        let response = match native.interaction {
            Interaction::RoughConductor {
                optics, surface, ..
            } => incoming.evaluate_smith_conductor_connection(
                cx,
                contact.witness.normal,
                &camera.scene.conductors[optics],
                surface,
                source_frame,
                bundle.interface_state(),
                settings,
                samples,
            ),
            Interaction::BrushedConductor { optics, detail, .. } => {
                let instance = camera
                    .scene
                    .geometry
                    .instances()
                    .get(contact.witness.boundary_id)
                    .ok_or(AreaLightCameraError::PathMismatch { binding })?;
                let transform = instance.transform();
                let local = transform
                    .inverse()
                    .map_err(SceneError::from)
                    .map_err(CurvedSceneError::from)
                    .map_err(CameraRadianceError::from)?
                    .transform_point(contact.witness.surface_point);
                incoming.evaluate_instanced_smith_brushed_conductor_connection(
                    cx,
                    &camera.scene.conductors[optics],
                    &detail,
                    transform,
                    local,
                    contact.witness.normal,
                    source_frame,
                    bundle.interface_state(),
                    settings,
                    samples,
                )
            }
            _ => unreachable!("continuous conductor admission above"),
        }
        .map_err(smith_error)?;
        Ok((
            face_normal(contact, endpoint.ray.dir),
            response.factor_per_sr(),
        ))
    }
}

fn smith_error(error: SmithBeamError) -> AreaLightCameraError {
    CameraRadianceError::from(CurvedSceneError::from(SceneError::from(error))).into()
}

//! Shared native receiving-vertex response for clear and participating paths.
//!
//! The source connection ends before the receiver scatters. Its native contact
//! determines the physical footprint and local material frame; the target's
//! shooting section is never used as a shading normal. No continuation draw or
//! independently sampled microfacet is introduced here.

use super::{AreaLightCameraError, diffuse_contribution, face_normal};
use super::super::super::{CameraRadianceError, CurvedCameraScene, phase_index};
use crate::dielectric::DielectricError;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::{Interaction, SceneError, bind_dielectric_surface};
use crate::polarization::dielectric::scene::curved::{CurvedContact, CurvedPathResult, CurvedSceneError};
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleConnection, BundleConnectionError, BundleConnectionTarget,
};
use crate::polarization::dielectric::transport::{PolarizationFrame, PolarizedBeam};
use crate::polarization::dielectric::transport::connection::SurfaceConnectionResponse;
use fs_geom::{Point3, Vec3};
use fs_math::det;

#[cfg(test)]
mod anisotropic_tests;

/// The two original adjacent media admitted by a continuous native receiver.
/// Fixed-medium surfaces repeat the same index; smooth delta laws have no
/// continuous solid-angle proposal and cannot become direct-light receivers.
pub(in super::super::super) fn receiver_media(interaction: Interaction) -> Option<[usize; 2]> {
    match interaction {
        Interaction::DepolarizingDiffuse { medium, .. }
        | Interaction::BrushedConductor { medium, .. } => Some([medium, medium]),
        Interaction::RoughConductor { medium, surface, .. } if !surface.is_delta() =>
            Some([medium, medium]),
        Interaction::Dielectric { outside, inside, surface } if !surface.is_delta() =>
            Some([outside, inside]),
        _ => None,
    }
}

/// Copy-only anchor state. Keeping the arriving beam preserves the native
/// polarization analyser and its local medium before rough transmission.
#[derive(Clone, Copy)]
pub(in super::super::super) struct SurfaceReceiver {
    pub(in super::super) contact: CurvedContact,
    incoming: PolarizedBeam,
    medium: usize,
    diffuse: Option<(Vec3, f64)>,
}

#[derive(Clone, Copy)]
pub(in super::super::super) struct SurfaceResponse {
    pub(in super::super::super) normal: Vec3,
    /// Diffuse reflectance for the legacy diffuse path; otherwise a full
    /// reciprocal polarized radiance-BSDF contraction per steradian.
    pub(in super::super::super) factor: f64,
    pub(in super::super) pdf: f64,
}

impl SurfaceReceiver {
    pub(in super::super::super) fn at_arrival(
        camera: &CurvedCameraScene<'_, '_>,
        path: &CurvedPathResult,
    ) -> Result<Option<Self>, AreaLightCameraError> {
        let contact = *path.last_contact().ok_or(AreaLightCameraError::NumericalRange)?;
        Self::at_contact(camera, contact, *path.beam(), path.medium_index())
    }

    pub(in super::super) fn at_contact(
        camera: &CurvedCameraScene<'_, '_>,
        contact: CurvedContact,
        incoming: PolarizedBeam,
        medium: usize,
    ) -> Result<Option<Self>, AreaLightCameraError> {
        let binding = camera.scene.bindings.get(contact.witness.boundary_id)
            .ok_or(AreaLightCameraError::NumericalRange)?;
        if binding.object_id != contact.object_id {
            return Err(AreaLightCameraError::NumericalRange);
        }
        if receiver_media(binding.interaction).is_none() { return Ok(None); }
        if arrival_medium(binding.interaction, contact, incoming.ray().dir)? != medium {
            return Err(AreaLightCameraError::NumericalRange);
        }
        let diffuse = if let Interaction::DepolarizingDiffuse { reflectance, .. } = binding.interaction {
            let rho = camera.scene.reflectances[reflectance].at(incoming.wavelength_nm())
                .map_err(|_| CameraRadianceError::Optics(DielectricError::InvalidWavelength))?;
            if rho == 0.0 { return Ok(None); }
            Some((face_normal(contact, incoming.ray().dir), rho))
        } else { None };
        Ok(Some(Self { contact, incoming, medium, diffuse }))
    }

    pub(in super::super::super) fn is_rough(self) -> bool { self.diffuse.is_none() }

    pub(in super::super::super) fn target(self, tolerance: f64)
        -> Result<BundleConnectionTarget, BundleConnectionError> {
        if let Some((normal, _)) = self.diffuse {
            BundleConnectionTarget::try_diffuse_receiver(self.contact.object_id,
                self.contact.witness.surface_point, normal, tolerance)
        } else {
            BundleConnectionTarget::try_surface_receiver(self.contact.object_id,
                self.contact.witness.surface_point,
                face_normal(self.contact, self.incoming.ray().dir), tolerance)
        }
    }

    pub(in super::super::super) fn response(
        self, camera: &CurvedCameraScene<'_, '_>, connection: &BundleConnection, binding: usize,
    ) -> Result<SurfaceResponse, AreaLightCameraError> {
        let path = connection.path();
        let bundle = path.bundle();
        let endpoint = bundle.propagation().propagation();
        if let Some((normal, rho)) = self.diffuse {
            let incidence = -normal.dot(endpoint.ray.dir);
            if !incidence.is_finite() { return Err(AreaLightCameraError::NumericalRange); }
            return Ok(SurfaceResponse { normal,
                factor: if incidence > 0.0 { rho } else { 0.0 },
                pdf: incidence.max(0.0) / core::f64::consts::PI });
        }
        let contact = *path.last_contact().ok_or(AreaLightCameraError::PathMismatch { binding })?;
        if !path.is_surface_receiver() || !same_contact(self.contact, contact) {
            return Err(AreaLightCameraError::PathMismatch { binding });
        }
        let frame = PolarizationFrame::try_new(endpoint.ray.dir, endpoint.transverse_axis)
            .map_err(CurvedSceneError::from).map_err(CameraRadianceError::from)?;
        let response = self.evaluate(camera, contact, endpoint.ray.origin, path.medium_index(),
            frame, bundle.interface_state(), binding)?;
        Ok(SurfaceResponse { normal: face_normal(contact, endpoint.ray.dir),
            factor: response.factor, pdf: response.pdf })
    }

    /// Validate the native recorded proposal at its exact sampled direction.
    /// A narrowly concentrated GGX lobe has no cosine-law Lipschitz bound, so
    /// a nearby numerically reconstructed source direction cannot replace it.
    pub(in super::super) fn sampled_pdf(
        self, camera: &CurvedCameraScene<'_, '_>, outgoing: Vec3,
        outgoing_medium: usize, binding: usize,
    ) -> Result<f64, AreaLightCameraError> {
        if let Some((normal, _)) = self.diffuse {
            return Ok(normal.dot(outgoing).max(0.0) / core::f64::consts::PI);
        }
        let arrival = outgoing.scale(-1.0);
        let axis = stable_axis(arrival);
        let frame = PolarizationFrame::try_new(arrival, axis)
            .map_err(CurvedSceneError::from).map_err(CameraRadianceError::from)?;
        Ok(self.evaluate(camera, self.contact, self.incoming.ray().origin, outgoing_medium,
            frame, Stokes::unpolarized(1.0), binding)?.pdf)
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate(
        self, camera: &CurvedCameraScene<'_, '_>, contact: CurvedContact, position: Point3,
        source_medium: usize, source_frame: PolarizationFrame, source_stokes: Stokes,
        binding: usize,
    ) -> Result<SurfaceConnectionResponse, AreaLightCameraError> {
        let native = camera.scene.bindings.get(contact.witness.boundary_id)
            .ok_or(AreaLightCameraError::PathMismatch { binding })?;
        if native.object_id != contact.object_id
            || arrival_medium(native.interaction, contact, self.incoming.ray().dir)? != self.medium
            || arrival_medium(native.interaction, contact, source_frame.direction())? != source_medium {
            return Err(AreaLightCameraError::PathMismatch { binding });
        }
        let wavelength = self.incoming.wavelength_nm();
        // Reconstruct only the local evaluation point. The original analyser,
        // direction and physical prefix are retained without extra propagation.
        let incoming = camera.scene.launch_beam(position, self.incoming.frame(), wavelength,
            self.medium, self.incoming.power()).map_err(CameraRadianceError::from)?;
        match native.interaction {
            Interaction::Dielectric { outside, inside, surface } if !surface.is_delta() => {
                let opposite = if self.medium == outside { inside } else { outside };
                let medium = camera.scene.medium_at(opposite, position, wavelength)
                    .map_err(CameraRadianceError::from)?;
                let instance = camera.scene.geometry.instances().get(contact.witness.boundary_id)
                    .ok_or(AreaLightCameraError::PathMismatch { binding })?;
                let surface = bind_dielectric_surface(surface, instance.transform())
                    .map_err(CurvedSceneError::from).map_err(CameraRadianceError::from)?;
                let response = if let Some(stack) = camera.scene.coating_stacks
                    .get(contact.witness.boundary_id).copied().flatten() {
                    incoming.evaluate_stack_dielectric_connection(contact.witness.normal, medium,
                        surface, stack, source_frame, source_stokes)
                } else {
                    incoming.evaluate_dielectric_connection(contact.witness.normal, medium, surface,
                        source_frame, source_stokes)
                };
                response.map_err(CurvedSceneError::from).map_err(CameraRadianceError::from).map_err(Into::into)
            }
            Interaction::RoughConductor { optics, surface, .. } if !surface.is_delta() =>
                incoming.evaluate_rough_conductor_connection(contact.witness.normal,
                    &camera.scene.conductors[optics], surface, source_frame, source_stokes)
                    .map_err(SceneError::from).map_err(CurvedSceneError::from)
                    .map_err(CameraRadianceError::from).map_err(Into::into),
            Interaction::BrushedConductor { optics, detail, .. } => {
                let instance = camera.scene.geometry.instances().get(contact.witness.boundary_id)
                    .ok_or(AreaLightCameraError::PathMismatch { binding })?;
                let transform = instance.transform();
                let local = transform.inverse().map_err(SceneError::from)
                    .map_err(CurvedSceneError::from).map_err(CameraRadianceError::from)?
                    .transform_point(contact.witness.surface_point);
                incoming.evaluate_instanced_brushed_conductor_connection(
                    &camera.scene.conductors[optics], &detail, transform, local,
                    contact.witness.normal, source_frame, source_stokes)
                    .map_err(SceneError::from).map_err(CurvedSceneError::from)
                    .map_err(CameraRadianceError::from).map_err(Into::into)
            }
            _ => Err(AreaLightCameraError::PathMismatch { binding }),
        }
    }

    pub(in super::super) fn contribution(
        self, initial_index: f64, factor: f64, irradiance: f64, weight: f64,
    ) -> Result<f64, AreaLightCameraError> {
        let local_index = phase_index(self.incoming.medium(), self.incoming.wavelength_nm())?;
        if self.diffuse.is_some() {
            diffuse_contribution(self.incoming.power().i, initial_index, local_index,
                factor, irradiance, weight)
        } else {
            directional_contribution(self.incoming.power().i, initial_index, local_index,
                factor, irradiance, weight)
        }
    }
}

pub(in super::super) fn directional_contribution(
    prefix: f64, initial_index: f64, local_index: f64,
    factor: f64, irradiance: f64, weight: f64,
) -> Result<f64, AreaLightCameraError> {
    let ratio = initial_index / local_index;
    let value = prefix * (ratio * ratio) * (factor * irradiance) * weight;
    if value.is_finite() && value > 0.0 { return Ok(value); }
    let value = det::exp(det::ln(prefix) + 2.0 * (det::ln(initial_index) - det::ln(local_index))
        + det::ln(factor) + det::ln(irradiance) + det::ln(weight));
    if value.is_finite() && value > 0.0 { Ok(value) } else { Err(AreaLightCameraError::NumericalRange) }
}

fn arrival_medium(interaction: Interaction, contact: CurvedContact, direction: Vec3)
    -> Result<usize, AreaLightCameraError> {
    match interaction {
        Interaction::DepolarizingDiffuse { medium, .. }
        | Interaction::RoughConductor { medium, .. }
        | Interaction::BrushedConductor { medium, .. } => Ok(medium),
        Interaction::Dielectric { outside, inside, .. } => {
            let cosine = contact.witness.normal.dot(direction);
            if !cosine.is_finite() || cosine == 0.0 {
                return Err(AreaLightCameraError::NumericalRange);
            }
            Ok(if cosine < 0.0 { outside } else { inside })
        }
        _ => Err(AreaLightCameraError::NumericalRange),
    }
}

fn same_contact(a: CurvedContact, b: CurvedContact) -> bool {
    a.object_id == b.object_id && a.geometry_identity == b.geometry_identity
        && a.frame_identity == b.frame_identity
        && a.witness.boundary_id == b.witness.boundary_id
        && a.witness.feature_id == b.witness.feature_id
}

fn stable_axis(direction: Vec3) -> Vec3 {
    let [x, y, z] = [direction.x.abs(), direction.y.abs(), direction.z.abs()];
    if x <= y && x <= z { Vec3::new(1.0, 0.0, 0.0) }
    else if y <= z { Vec3::new(0.0, 1.0, 0.0) }
    else { Vec3::new(0.0, 0.0, 1.0) }
}

//! Prepared scalar wave imaging in the original clear homogeneous medium.
//!
//! The original camera owns pupil diffraction, material-relative lens power
//! and physical medium wavelength. This owner proves that the original opaque
//! first hits supply its geometric axial source depths: all forward meshes
//! lie on the declared planes, and the existing full-pupil footprint proof
//! excludes inter-plane clipping. Every spectral lane still traces the native
//! source scene and retains its own source spectrum, contact and optical path.
//!
//! This is the camera's existing scalar paraxial incoherent model, including
//! its stated finite phase and angular limits. It does not infer object depth
//! from a curved optical path, model coherent scene interference, or account
//! for obstacle edge diffraction. Refracting/reflecting primary surfaces,
//! gradient-index propagation, volume emission, environments and absorbing
//! pupil apodization require a different source/pupil model and are refused.

use super::sensor::CurvedSensorPixelSample;
use super::*;
use crate::camera::{LayeredDiffraction, ProjectedPupilFootprint};
use crate::dielectric::BeerLambertParameters;
use crate::instances::SharedGeometry;
use crate::polarization::dielectric::scene::curved::SceneMedium;
use crate::polarization::dielectric::scene::{Interaction, PathTermination};
use crate::spectral::detector::QuantumEfficiencySensor;
use fs_geom::Point3;

/// A failed proof leaves no prepared wave observer or partial pixel result.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WaveSceneError {
    InvalidMedium,
    AbsorbingMedium,
    ExteriorMismatch,
    UnsupportedSource,
    UnsupportedSurface { object_id: u64 },
    UnsupportedGeometry { object_id: u64 },
    SourcePlane { object_id: u64 },
    PupilOcclusion { nearer: usize, farther: usize },
    UnresolvedPrimary,
}

/// An immutable camera/source binding for clear homogeneous wave imaging.
///
/// Construction proves geometric full-pupil visibility over the original
/// rigid mesh instances. The original medium law must equal the camera's
/// explicit exterior law; an implicit vacuum camera accepts only a genuinely
/// nondispersive index-one law. Zero extinction must be explicitly declared.
/// Wavelength support, phase/depth admission and emission are checked at every
/// requested native lane before any XYZ or sensor result is returned.
///
/// The private raster cannot be reused with another source scene or medium.
/// Samples allocate no observer-owned workspace and retain the existing
/// independent raster/lens/spectral/event counter streams and path budgets.
pub struct CurvedWaveCameraRaster<'a, 's, 'm> {
    raster: CurvedCameraRaster<'a>,
    scene: &'a CurvedCameraScene<'s, 'm>,
    initial_medium: usize,
}

impl<'a, 's, 'm> CurvedWaveCameraRaster<'a, 's, 'm> {
    /// Bind finite single-plane defocus or layered diffraction to the original
    /// opaque source geometry and explicitly selected homogeneous medium.
    /// The whole virtual pupil lies in that same clear medium; no enclosure
    /// or equality at one reference wavelength is used to infer this fact.
    pub fn try_new(
        cx: &Cx<'_>,
        camera: &'a PhysicalCamera,
        scene: &'a CurvedCameraScene<'s, 'm>,
        initial_medium: usize,
        width: u32,
        height: u32,
    ) -> Result<Self, CurvedRasterError> {
        cx.checkpoint().map_err(CameraError::from)?;
        let raster = CurvedCameraRaster::prepare(camera, width, height)?;
        admit_medium(camera, scene, initial_medium)?;
        admit_sources(cx, camera, scene)?;
        Ok(Self {
            raster,
            scene,
            initial_medium,
        })
    }

    pub const fn width(&self) -> u32 {
        self.raster.width()
    }
    pub const fn height(&self) -> u32 {
        self.raster.height()
    }

    /// Observe one absolute pixel/sample identity, applying the original
    /// spectral pupil weights exactly once before visible XYZ integration.
    pub fn sample_pixel(
        &self,
        cx: &Cx<'_>,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedPixelSample, CurvedRasterError> {
        self.raster.sample_pixel(
            cx,
            self.scene,
            self.initial_medium,
            settings,
            seed,
            pixel,
            sample,
        )
    }

    /// Observe the same admitted wave image through the original sum-of-QE
    /// proposal and channel integrals. Photon energy uses vacuum wavelength;
    /// the medium wavelength belongs only to the camera's pupil propagation.
    pub fn sample_pixel_sensor<const C: usize>(
        &self,
        cx: &Cx<'_>,
        sensor: &QuantumEfficiencySensor<C>,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedSensorPixelSample<C>, CurvedRasterError> {
        self.raster.sample_pixel_sensor(
            cx,
            self.scene,
            sensor,
            self.initial_medium,
            settings,
            seed,
            pixel,
            sample,
        )
    }
}

fn admit_medium(
    camera: &PhysicalCamera,
    source: &CurvedCameraScene<'_, '_>,
    medium: usize,
) -> Result<(), CurvedRasterError> {
    let Some(SceneMedium::Homogeneous(material)) = source.scene.media.get(medium) else {
        return Err(WaveSceneError::InvalidMedium.into());
    };
    if !matches!(
        material.absorption().parameters(),
        BeerLambertParameters::Clear
    ) {
        return Err(WaveSceneError::AbsorbingMedium.into());
    }
    let ior = material.ior();
    if let Some(exterior) = camera
        .dispersive_lens()
        .and_then(|lens| lens.exterior_ior())
    {
        if exterior != ior {
            return Err(WaveSceneError::ExteriorMismatch.into());
        }
    } else if ior.is_dispersive() || ior.eval(ior.wavelength_bounds_nm().0).ok() != Some(1.0) {
        return Err(WaveSceneError::ExteriorMismatch.into());
    }
    Ok(())
}

fn admit_sources(
    cx: &Cx<'_>,
    camera: &PhysicalCamera,
    source: &CurvedCameraScene<'_, '_>,
) -> Result<(), CurvedRasterError> {
    if source.environment.is_some()
        || source.volume_emitters.iter().any(Option::is_some)
        || source.spatial_volume_emitters.iter().any(Option::is_some)
    {
        return Err(WaveSceneError::UnsupportedSource.into());
    }
    let layered = camera.layered_diffraction();
    let single = camera
        .defocused_diffraction()
        .and_then(|p| p.object_distance_m());
    let single_depths = [single.unwrap_or(0.0)];
    let depths = if let Some(policy) = &layered {
        policy.object_distances_m()
    } else if single.is_some() {
        &single_depths
    } else {
        return Err(CameraError::InvalidDiffraction.into());
    };
    let mut planes: [Option<ProjectedPupilFootprint>; LayeredDiffraction::MAX_PLANES] =
        [None; LayeredDiffraction::MAX_PLANES];
    let eye = camera.eye();
    for instance in source.scene.geometry.instances() {
        cx.checkpoint().map_err(CameraError::from)?;
        let object_id = instance.object_id();
        let binding = source
            .scene
            .bindings
            .binary_search_by_key(&object_id, |b| b.object_id)
            .map_err(|_| WaveSceneError::UnsupportedSurface { object_id })?;
        if source.scene.bindings[binding].interaction != Interaction::Opaque {
            return Err(WaveSceneError::UnsupportedSurface { object_id }.into());
        }
        let SharedGeometry::Mesh(mesh) = instance.geometry() else {
            return Err(WaveSceneError::UnsupportedGeometry { object_id }.into());
        };
        let transform = instance.transform();
        let origin = transform.transform_point(Point3::new(0.0, 0.0, 0.0));
        let mut behind_pupil = None;
        let mut primitive_plane = None;
        for (vertex_index, &vertex) in mesh.vertices.iter().enumerate() {
            if vertex_index % 256 == 0 {
                cx.checkpoint().map_err(CameraError::from)?;
            }
            let point = transform.transform_point(Point3::new(vertex[0], vertex[1], vertex[2]));
            let displacement = point.delta_from(eye);
            let depth = displacement.dot(camera.forward());
            // Include original local coordinates and placement as well as the
            // resulting world point: a large cancelling translation must not
            // masquerade as a precisely represented small source distance.
            let scale = [
                point.x, point.y, point.z, eye.x, eye.y, eye.z, origin.x, origin.y, origin.z,
                vertex[0], vertex[1], vertex[2],
            ]
            .into_iter()
            .map(f64::abs)
            .fold(1.0_f64, f64::max);
            let roundoff = 128.0 * f64::EPSILON * scale;
            if !depth.is_finite() || !roundoff.is_finite() || depth.abs() <= roundoff {
                return Err(WaveSceneError::SourcePlane { object_id }.into());
            }
            let behind = depth < 0.0;
            if behind_pupil.is_some_and(|previous| previous != behind) {
                return Err(WaveSceneError::SourcePlane { object_id }.into());
            }
            behind_pupil = Some(behind);
            if behind {
                continue;
            }
            let mut matched = None;
            for (layer, &declared) in depths.iter().enumerate() {
                let tolerance = roundoff + 128.0 * f64::EPSILON * declared.max(1.0);
                if (depth - declared).abs() <= tolerance {
                    if matched.is_some() {
                        return Err(WaveSceneError::SourcePlane { object_id }.into());
                    }
                    matched = Some(layer);
                }
            }
            let layer = matched.ok_or(WaveSceneError::SourcePlane { object_id })?;
            if primitive_plane.is_some_and(|previous| previous != layer) {
                return Err(WaveSceneError::SourcePlane { object_id }.into());
            }
            primitive_plane = Some(layer);
            let footprint = ProjectedPupilFootprint::point(
                [
                    displacement.dot(camera.right()),
                    displacement.dot(camera.up()),
                    depth,
                ],
                roundoff,
            )
            .ok_or(WaveSceneError::SourcePlane { object_id })?;
            if let Some(existing) = &mut planes[layer] {
                existing.include(footprint);
            } else {
                planes[layer] = Some(footprint);
            }
        }
    }
    let radius = camera.aperture().radius_m();
    for nearer in 0..depths.len() {
        cx.checkpoint().map_err(CameraError::from)?;
        let Some(front) = planes[nearer] else {
            continue;
        };
        for (farther, back) in planes
            .iter()
            .enumerate()
            .take(depths.len())
            .skip(nearer + 1)
        {
            if back.is_some_and(|back| {
                !front.clears(back, [(-radius).next_down(); 2], [radius.next_up(); 2])
            }) {
                return Err(WaveSceneError::PupilOcclusion { nearer, farther }.into());
            }
        }
    }
    cx.checkpoint().map_err(CameraError::from)?;
    Ok(())
}

// Opaque first hits retain the original geometric witness. Its source point
// supplies axial object depth; the integrated beam endpoint and n*ds/group
// optical lengths remain untouched in the retained native path.
pub(super) fn source_point(
    sample: &CameraRadianceSample,
) -> Result<Option<Point3>, CurvedRasterError> {
    let path = sample.path();
    match sample.endpoint() {
        super::super::RadianceEndpoint::Surface { object_id }
            if path.termination() == PathTermination::Opaque && path.surface_events() == 1 =>
        {
            let contact = path
                .last_contact()
                .ok_or(WaveSceneError::UnresolvedPrimary)?;
            if contact.object_id != object_id {
                return Err(WaveSceneError::UnresolvedPrimary.into());
            }
            Ok(Some(contact.witness.surface_point))
        }
        super::super::RadianceEndpoint::Environment
            if path.termination() == PathTermination::Escaped
                && path.surface_events() == 0
                && sample.radiance_w_m2_sr_nm() == 0.0 =>
        {
            Ok(None)
        }
        _ => Err(WaveSceneError::UnresolvedPrimary.into()),
    }
}

#[cfg(test)]
mod tests;

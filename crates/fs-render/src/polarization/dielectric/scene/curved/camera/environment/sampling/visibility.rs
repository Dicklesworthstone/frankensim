//! First-boundary visibility in a declared straight, nonabsorbing medium.
//!
//! The real receiving contact owns surface departure; volume vertices have
//! no artificial contact or positional offset. A finite integration budget
//! does not stand in for infinity: the original outward-bounded scene oracle
//! must also prove that the remaining forward halfline misses every surface.

use super::super::super::CameraRadianceError;
use super::{EnvironmentCameraError, EnvironmentCameraSettings};
use crate::charts::Ray;
use crate::gradient_index::HomogeneousIndexField;
use crate::gradient_index::event::{
    BoundarySide, BoundaryTermination, GradientDeparture, propagate_to_boundary_with_control,
};
use crate::polarization::dielectric::scene::curved::{
    CurvedContact, CurvedPolarizedScene, CurvedSceneError, SceneMedium,
};
use fs_exec::Cx;
use fs_geom::Vec3;

#[derive(Clone, Copy, Default)]
pub(super) struct Visibility {
    pub(super) visible: bool,
    pub(super) fields: u64,
    pub(super) geometry: u64,
}

pub(super) fn supports_segment(
    scene: &CurvedPolarizedScene<'_>,
    medium: usize,
    wavelength: f64,
) -> Result<bool, EnvironmentCameraError> {
    match scene
        .media
        .get(medium)
        .ok_or(EnvironmentCameraError::NumericalRange)?
    {
        SceneMedium::Homogeneous(material) => Ok(material
            .absorption()
            .extinction_per_m(wavelength)
            .map_err(CameraRadianceError::from)?
            == 0.0),
        // A sampled zero gradient does not prove a straight angular map.
        // These paths retain their original camera estimator without NEE.
        SceneMedium::GradientIndex { .. } | SceneMedium::Differentiable { .. } => Ok(false),
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn clear_to_environment(
    cx: &Cx<'_>,
    scene: &CurvedPolarizedScene<'_>,
    contact: Option<CurvedContact>,
    ray: Ray,
    transverse_axis: Vec3,
    wavelength: f64,
    medium: usize,
    settings: EnvironmentCameraSettings,
    remaining_work: u64,
) -> Result<Visibility, EnvironmentCameraError> {
    cx.checkpoint()?;
    if !supports_segment(scene, medium, wavelength)? {
        return Err(EnvironmentCameraError::NumericalRange);
    }
    let SceneMedium::Homogeneous(material) = scene.media[medium] else {
        return Err(EnvironmentCameraError::NumericalRange);
    };
    if let Some(contact) = contact {
        let instance = scene
            .geometry
            .instances()
            .get(contact.witness.boundary_id)
            .ok_or(EnvironmentCameraError::NumericalRange)?;
        if instance.object_id() != contact.object_id
            || instance.geometry_identity() != contact.geometry_identity
            || instance.frame_identity() != contact.frame_identity
        {
            return Err(EnvironmentCameraError::NumericalRange);
        }
    }
    let tolerance = settings.visibility.surface_tolerance_m;
    if scene
        .boundaries
        .straight_ray_escapes(&ray, tolerance)
        .map_err(CurvedSceneError::from)
        .map_err(CameraRadianceError::from)?
    {
        return Ok(Visibility {
            visible: true,
            ..Visibility::default()
        });
    }
    if remaining_work == 0 {
        return Err(EnvironmentCameraError::GeometryLimit);
    }
    let departure = if let Some(contact) = contact {
        let cosine = contact.witness.normal.dot(ray.dir);
        if !cosine.is_finite() || cosine == 0.0 {
            return Err(EnvironmentCameraError::NumericalRange);
        }
        Some(GradientDeparture {
            witness: contact.witness,
            side: if cosine > 0.0 {
                BoundarySide::Positive
            } else {
                BoundarySide::Negative
            },
        })
    } else {
        None
    };
    let mut oracle = scene.boundaries.oracle(cx, remaining_work, tolerance);
    let field = HomogeneousIndexField::new(material);
    let segment = propagate_to_boundary_with_control(
        &field,
        &mut oracle,
        ray,
        transverse_axis,
        wavelength,
        settings.visibility_distance_m,
        departure,
        settings.visibility,
        || cx.checkpoint().is_ok(),
    )
    .map_err(|error| {
        CameraRadianceError::from(oracle.last_error().map_or(
            CurvedSceneError::Propagation(error),
            CurvedSceneError::Geometry,
        ))
    })?;
    let visible = match segment.termination {
        BoundaryTermination::Surface(_) => false,
        BoundaryTermination::DistanceLimit => {
            if !scene
                .boundaries
                .straight_ray_escapes(&segment.propagation.ray, tolerance)
                .map_err(CurvedSceneError::from)
                .map_err(CameraRadianceError::from)?
            {
                return Err(EnvironmentCameraError::UnresolvedVisibility);
            }
            true
        }
    };
    cx.checkpoint()?;
    Ok(Visibility {
        visible,
        fields: u64::from(segment.propagation.field_evaluations),
        geometry: oracle.work_used(),
    })
}

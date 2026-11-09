//! Fixed-pupil material-dispersive motion on the original correspondence.
use super::*;
use crate::camera::LensSample;

/// Project the same admitted material point at one actual geometric-lens
/// wavelength and pupil sample across all three frames. Identity, cuts,
/// normals, raster signs and target visibility remain with their original
/// owners. Each distinct projection uses one original ray and no iteration.
#[allow(clippy::too_many_arguments)]
pub fn compute_geometric_dispersive_motion_vectors(
    cx: &Cx<'_>,
    sample: PrimarySurfaceSample,
    previous: &MotionFrame,
    current: &MotionFrame,
    next: &MotionFrame,
    raster: RasterSize,
    wavelength_nm: f64,
    lens_sample: LensSample,
    position_tolerance_m: f64,
) -> Result<MotionVectorComputation, MotionVectorError> {
    cx.checkpoint()?;
    let mut cached = None;
    let result = compute_motion_vectors_with(
        sample,
        previous,
        current,
        next,
        raster,
        |frame, point, raster| {
            let is_current = core::ptr::eq(frame, current);
            if is_current {
                if let Some(projection) = cached {
                    return Ok(projection);
                }
            }
            let camera = frame.camera();
            let depth = camera.forward().dot(point.delta_from(camera.eye()));
            if !depth.is_finite() {
                return Err(MotionVectorError::InvalidProjection);
            }
            if depth <= 0.0 {
                return Ok(ProjectedWorldPoint::BehindCamera {
                    signed_depth_m: depth,
                });
            }
            let native = camera
                .project_geometric_spectral_point(
                    cx,
                    point,
                    raster.aspect_ratio(),
                    wavelength_nm,
                    lens_sample,
                    position_tolerance_m,
                )
                .map_err(|error| match error {
                    CameraError::Cancelled => MotionVectorError::Cancelled,
                    other => MotionVectorError::Camera(other),
                })?;
            let projection = raster_projection(native.ndc_xy(), native.depth_m(), raster)?;
            if is_current {
                cached = Some(projection);
            }
            Ok(projection)
        },
    )?;
    cx.checkpoint()?;
    Ok(result)
}

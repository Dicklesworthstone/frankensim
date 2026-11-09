//! Conditional geometric motion through the original sequential lens.
use super::*;
use crate::camera::{LensSample, SequentialProjectionError, SequentialProjectionOptions};

impl From<SequentialProjectionError> for MotionVectorError {
    fn from(error: SequentialProjectionError) -> Self {
        match error {
            SequentialProjectionError::Camera(CameraError::Cancelled) => Self::Cancelled,
            _ => Self::SequentialProjection(error),
        }
    }
}

/// Follow one stable material point through a fixed wavelength and pupil draw.
///
/// The current seed selects a local lens-image branch. Its converged tangent
/// coordinate seeds both reference frames; each frame uses its original native
/// lens prescription, coatings and ghost selection. This conditional vector is
/// neither a pupil average nor a visibility decision. Use the existing target
/// observation/reprojection validation before reusing a target-frame sample.
///
/// All correspondence, identity, normal, camera-cut and raster semantics are
/// shared with `compute_motion_vectors`. The work limit applies separately to
/// each distinct projection. Passing `current` as an unused endpoint reuses its
/// accepted projection without tracing it again. A clipped or unresolved
/// inverse returns a typed error and never a fabricated zero vector.
#[allow(clippy::too_many_arguments)]
pub fn compute_sequential_motion_vectors(
    cx: &Cx<'_>,
    sample: PrimarySurfaceSample,
    previous: &MotionFrame,
    current: &MotionFrame,
    next: &MotionFrame,
    raster: RasterSize,
    wavelength_nm: f64,
    lens_sample: LensSample,
    mut options: SequentialProjectionOptions,
) -> Result<MotionVectorComputation, MotionVectorError> {
    cx.checkpoint()?;
    let mut current_projection = None;
    let result = compute_motion_vectors_with(
        sample,
        previous,
        current,
        next,
        raster,
        |frame, point, raster| {
            let is_current = core::ptr::eq(frame, current);
            if is_current {
                if let Some(projection) = current_projection {
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
            let native = camera.project_sequential_point(
                cx,
                point,
                raster.aspect_ratio(),
                wavelength_nm,
                lens_sample,
                options,
            )?;
            let projection = raster_projection(native.ndc_xy(), native.depth_m(), raster)?;
            if is_current {
                options.initial_tangent_xy = native.tangent_xy();
                current_projection = Some(projection);
            }
            Ok(projection)
        },
    )?;
    cx.checkpoint()?;
    Ok(result)
}

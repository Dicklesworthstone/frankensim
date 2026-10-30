//! Accepted surface, source and geometry states over one physical exposure.
//!
//! Frames borrow complete accepted scenes. Their lighting proposals are prepared
//! once, then shared by all paths selecting that state. Temperature, emissivity,
//! optical materials and topology are held over their declared intervals; no
//! interpolation invents intermediate simulation states or changes owner energy.

use super::*;
#[cfg(feature = "volumes")]
use crate::volumes::transport::timeline::MediumTimeline;

mod execution;
pub use execution::{
    SceneExposureFilm, render_scene_exposure, render_scene_exposure_range,
    render_scene_exposure_range_with_execution, render_scene_exposure_with_execution,
};

/// One immutable accepted scene held over a finite, nonzero physical interval.
///
/// Intervals are `[start_s, end_s)`: the later state wins at a discontinuity.
/// Only the final history endpoint is inclusive. A frame may contain its own
/// continuous rigid trajectories, evaluated at the same absolute physical time.
#[derive(Clone, Copy)]
pub struct SceneFrame<'a> {
    /// Inclusive accepted state start, in absolute physical seconds.
    pub start_s: f64,
    /// Exclusive accepted state end, except at the final history endpoint.
    pub end_s: f64,
    /// Accepted geometry, materials, sources and scene camera, without copying.
    pub scene: &'a Scene,
}

/// Explicit bounds on the scene history and prepared light-source metadata.
///
/// Totals count all supplied frames, including repeated references to one scene.
/// Existing geometry and environment storage remains borrowed. Mesh continuity
/// hashing reads its existing vertices and faces with cancellation checkpoints.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SceneExposureBudget {
    /// Maximum accepted intervals inspected and bound to the exposure.
    pub max_frames: usize,
    /// Maximum sum of scene primitive counts across the supplied history.
    pub max_total_primitives: usize,
    /// Maximum source triangle count across the supplied history: each static
    /// or moving positive-emissivity thermal mesh face counts once, and each
    /// rectangular light counts as its two matching geometry triangles.
    pub max_total_light_triangles: usize,
}

/// Missing, ambiguous or excessive accepted scene history is refused.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SceneExposureError {
    /// At least one accepted scene interval is required.
    Empty,
    /// History metadata exceeds the caller's limit or cannot be allocated.
    Budget,
    /// An interval has nonfinite, equal or reversed endpoints.
    InvalidInterval {
        /// Invalid interval index.
        index: usize,
    },
    /// Adjacent frames leave time without an accepted scene.
    Gap {
        /// Index of the later frame.
        index: usize,
    },
    /// Adjacent frames overlap or are out of order.
    Overlap {
        /// Index of the later frame.
        index: usize,
    },
    /// A requested time is nonfinite or outside the admitted exposure.
    TimeOutsideExposure,
    /// The complete shutter is not covered by the accepted scene history.
    ShutterOutsideHistory,
}

impl core::fmt::Display for SceneExposureError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "accepted scene exposure refused: {self:?}")
    }
}

impl core::error::Error for SceneExposureError {}

impl From<SceneExposureError> for TracerError {
    fn from(error: SceneExposureError) -> Self {
        Self::SceneExposure(error)
    }
}

/// Contiguous accepted scenes and their once-prepared source proposals.
///
/// Every camera sample selects one scene at its physical shutter time and keeps
/// that scene and its matching light proposals for the complete path. Optional
/// cinematic cameras and physical media use that same clock. This is a
/// quasi-static exposure model, without finite light travel time.
///
/// The continuity guard binds actual geometry, material, source and camera
/// values. Direct chart backends lack a general content-addressing contract:
/// their existing process-local allocation guard detects replacement but is
/// not an authoritative persisted content identity.
pub struct SceneExposure<'a> {
    pub(super) frames: &'a [SceneFrame<'a>],
    pub(super) lighting: Vec<Option<SceneLighting<'a>>>,
    pub(super) shutter: ShutterInterval,
    pub(super) camera: CameraPath<'a>,
    #[cfg(feature = "volumes")]
    pub(super) medium: Option<MediumInput<'a, 'a>>,
    history_identity: ContentHash,
    camera_identity: Option<ContentHash>,
}

impl<'a> SceneExposure<'a> {
    /// Admit exact contiguous state coverage and prepare only selected frames.
    ///
    /// A scene's trajectories need cover only the intersection of that frame
    /// with this exposure. An exact shutter endpoint at a state change selects
    /// the later frame, including a zero-duration shutter. Admission never
    /// clamps missing support or interpolates across accepted state changes.
    pub fn try_new(
        cx: &Cx<'_>,
        frames: &'a [SceneFrame<'a>],
        shutter: ShutterInterval,
        budget: SceneExposureBudget,
    ) -> Result<Self, TracerError> {
        cx.checkpoint()?;
        if frames.is_empty() {
            return Err(SceneExposureError::Empty.into());
        }
        if frames.len() > budget.max_frames {
            return Err(SceneExposureError::Budget.into());
        }

        // Refuse the full budget before constructing any light hierarchy.
        let mut primitive_count = 0_usize;
        let mut light_triangles = 0_usize;
        for (index, frame) in frames.iter().enumerate() {
            cx.checkpoint()?;
            if !frame.start_s.is_finite()
                || !frame.end_s.is_finite()
                || frame.start_s >= frame.end_s
            {
                return Err(SceneExposureError::InvalidInterval { index }.into());
            }
            if index > 0 {
                let previous_end = frames[index - 1].end_s;
                if frame.start_s > previous_end {
                    return Err(SceneExposureError::Gap { index }.into());
                }
                if frame.start_s < previous_end {
                    return Err(SceneExposureError::Overlap { index }.into());
                }
            }
            add_budgeted(
                &mut primitive_count,
                frame.scene.primitives.len(),
                budget.max_total_primitives,
            )?;
            let rectangular_triangles = frame
                .scene
                .lights
                .len()
                .checked_mul(2)
                .ok_or(SceneExposureError::Budget)?;
            add_budgeted(
                &mut light_triangles,
                rectangular_triangles,
                budget.max_total_light_triangles,
            )?;
            for primitive in &frame.scene.primitives {
                cx.checkpoint()?;
                add_budgeted(
                    &mut light_triangles,
                    thermal_triangle_count(primitive),
                    budget.max_total_light_triangles,
                )?;
            }
        }
        if shutter.open_s() < frames[0].start_s
            || shutter.close_s() > frames[frames.len() - 1].end_s
        {
            return Err(SceneExposureError::ShutterOutsideHistory.into());
        }

        let mut lighting = Vec::new();
        lighting
            .try_reserve_exact(frames.len())
            .map_err(|_| SceneExposureError::Budget)?;
        let mut hash = DomainHasher::new("org.frankensim.render.scene-exposure-history.v1");
        hash.update(&(frames.len() as u64).to_le_bytes());
        for (index, frame) in frames.iter().enumerate() {
            cx.checkpoint()?;
            update_continuity_f64(&mut hash, frame.start_s);
            update_continuity_f64(&mut hash, frame.end_s);
            hash.update(scene_continuity_fingerprint(frame.scene, cx)?.as_bytes());
            let intersects = frame.start_s <= shutter.close_s()
                && (frame.end_s > shutter.open_s()
                    || (index + 1 == frames.len() && frame.end_s == shutter.open_s()));
            let prepared = if intersects {
                // Preserve exact accepted endpoint bits; reconstructing via a
                // nominal frame time can round a boundary onto another state.
                let clipped = ShutterInterval::try_from_canonical_parts(
                    frame.start_s.max(shutter.open_s()),
                    frame.end_s.min(shutter.close_s()),
                    shutter.convention(),
                    shutter.distribution(),
                )
                .map_err(|_| TracerError::InvalidInput)?;
                // Dark accepted states are meaningful. A physical medium may
                // also be attached after this admission and supply the light.
                Some(validate_scene_with_medium_emission(
                    frame.scene,
                    cx,
                    Some(clipped),
                    true,
                )?)
            } else {
                None
            };
            lighting.push(prepared);
        }
        cx.checkpoint()?;
        Ok(Self {
            frames,
            lighting,
            shutter,
            camera: CameraPath::Legacy,
            #[cfg(feature = "volumes")]
            medium: None,
            history_identity: hash.finalize(),
            camera_identity: None,
        })
    }

    /// Evaluate an admitted thin-lens/keyframed camera on the scene's clock.
    /// Existing shot coverage and camera-cut selection rules remain binding.
    pub fn with_camera(
        mut self,
        cx: &Cx<'_>,
        camera: &'a AnimatedCamera,
        cut_side: CutSide,
    ) -> Result<Self, TracerError> {
        let exposure = camera.admit_shutter(cx, self.shutter, cut_side)?;
        self.camera_identity = Some(camera_continuity_fingerprint(cx, camera)?);
        self.camera = CameraPath::Cinematic { camera, exposure };
        Ok(self)
    }

    /// Select immutable physical medium snapshots at the same sampled time.
    /// Surface and medium state changes may have different interval boundaries.
    #[cfg(feature = "volumes")]
    pub fn with_medium_history(
        mut self,
        cx: &Cx<'_>,
        history: &'a MediumTimeline<'a, 'a>,
    ) -> Result<Self, TracerError> {
        history.admit_shutter(cx, self.shutter)?;
        self.medium = Some(MediumInput::Timeline(history));
        Ok(self)
    }

    /// Admitted absolute exposure and its temporal sampling policy.
    #[must_use]
    pub const fn shutter(&self) -> ShutterInterval {
        self.shutter
    }

    /// Number of supplied accepted scene states.
    #[must_use]
    pub const fn frame_count(&self) -> usize {
        self.frames.len()
    }

    /// Numerical scene-history continuity guard, with the chart limitation
    /// documented on [`Self`]. This is not a persisted authority-bearing digest.
    #[must_use]
    pub const fn history_identity(&self) -> ContentHash {
        self.history_identity
    }

    pub(super) fn binding_identity(&self) -> ContentHash {
        let mut hash = DomainHasher::new("org.frankensim.render.scene-exposure-binding.v1");
        hash.update(self.history_identity.as_bytes());
        match self.camera_identity {
            Some(identity) => {
                hash.update(&[1]);
                hash.update(identity.as_bytes());
            }
            None => hash.update(&[0]),
        }
        #[cfg(feature = "volumes")]
        match self.medium {
            Some(MediumInput::Timeline(history)) => {
                hash.update(&[1]);
                hash.update(history.content_identity().as_bytes());
            }
            Some(MediumInput::Static(PathMedium::Absorbing(medium))) => {
                hash.update(&[2]);
                hash.update(medium.content_identity().as_bytes());
            }
            Some(MediumInput::Static(PathMedium::Scattering(medium))) => {
                hash.update(&[3]);
                hash.update(medium.content_identity().as_bytes());
            }
            None => hash.update(&[0]),
        }
        hash.finalize()
    }

    pub(super) fn initial_scene(&self) -> &'a Scene {
        self.frames[self.frame_index(self.shutter.open_s())].scene
    }

    pub(super) fn at_time(
        &self,
        absolute_time_s: f64,
    ) -> Result<(&'a Scene, &SceneLighting<'a>), TracerError> {
        if !absolute_time_s.is_finite()
            || absolute_time_s < self.shutter.open_s()
            || absolute_time_s > self.shutter.close_s()
        {
            return Err(SceneExposureError::TimeOutsideExposure.into());
        }
        let index = self.frame_index(absolute_time_s);
        let lighting = self.lighting[index]
            .as_ref()
            .ok_or(TracerError::InvalidInput)?;
        Ok((self.frames[index].scene, lighting))
    }

    fn frame_index(&self, absolute_time_s: f64) -> usize {
        // Admission guarantees a nonempty history covering every query here.
        // upper_bound(start) gives the later state at an exact discontinuity.
        self.frames
            .partition_point(|frame| frame.start_s <= absolute_time_s)
            .saturating_sub(1)
    }
}

fn add_budgeted(total: &mut usize, amount: usize, limit: usize) -> Result<(), TracerError> {
    *total = total
        .checked_add(amount)
        .filter(|&value| value <= limit)
        .ok_or(SceneExposureError::Budget)?;
    Ok(())
}

fn thermal_triangle_count(primitive: &Primitive) -> usize {
    let Some(surface) = primitive.material.thermal_surface() else {
        return 0;
    };
    if !surface.is_emitting() {
        return 0;
    }
    match &primitive.shape {
        Shape::Mesh(mesh) => mesh.triangles.len(),
        Shape::Instance(instance) => match instance.geometry() {
            SharedGeometry::Mesh(mesh) => mesh.triangles.len(),
            SharedGeometry::Chart(_) => 0,
        },
        Shape::AnimatedInstance(instance) => match instance.geometry() {
            SharedGeometry::Mesh(mesh) => mesh.triangles.len(),
            SharedGeometry::Chart(_) => 0,
        },
        Shape::Chart(_) => 0,
    }
}

#[cfg(test)]
mod tests;

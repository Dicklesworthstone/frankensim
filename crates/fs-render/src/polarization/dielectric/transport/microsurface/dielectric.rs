//! Two-sided Smith/GGX dielectric microsurfaces with native polarized optics.
//!
//! Heitz et al. (2016), section 7 and figure 11, reverse the shared statistical
//! height after transmission. In depth d=-ln(C1(h)), the opposite side is
//! -ln(1-exp(-d)); resetting depth or reusing the same depth loses the joint
//! geometry and biases glass. Both sides retain the same authored slope law.
//! <https://doi.org/10.1145/2897824.2925943>
//!
//! Each facet uses the original smooth Fresnel/coating and frame owner.
//! Transmission changes the local medium; total internal reflection retains its
//! phase. Statistical heights add no metric distance, delay or bulk absorption.
//! Optical film absorption is local to each facet, not an invented volume path.

use super::super::{BeamError, BeamMedium, PolarizationFrame, PolarizedBeam};
use super::{SmithError, SmithEvent, SmithGgx, SmithWalkSettings, unit};
use crate::dielectric::{
    DielectricCoatingStack, DielectricError, DielectricEvent, DielectricInterface,
    DielectricSurface,
};
use fs_exec::Cx;
use fs_geom::Vec3;
use fs_math::det;

/// Fixed-direction full-glass radiance response with native polarized facets.
pub mod connection;

/// A failed sample is never substituted with black or a truncated trajectory.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SmithDielectricError {
    Geometry(SmithError),
    Beam(BeamError),
    FacetLimit { facet_events: u32 },
    Cancelled,
}
impl From<SmithError> for SmithDielectricError {
    fn from(e: SmithError) -> Self {
        Self::Geometry(e)
    }
}
impl From<BeamError> for SmithDielectricError {
    fn from(e: BeamError) -> Self {
        Self::Beam(e)
    }
}
impl From<DielectricError> for SmithDielectricError {
    fn from(e: DielectricError) -> Self {
        Self::Beam(e.into())
    }
}
impl From<fs_exec::Cancelled> for SmithDielectricError {
    fn from(_: fs_exec::Cancelled) -> Self {
        Self::Cancelled
    }
}
impl core::fmt::Display for SmithDielectricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "polarized Smith dielectric: {self:?}")
    }
}
impl core::error::Error for SmithDielectricError {}

/// Completed ray-power path through either side of one statistical interface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmithDielectricSample {
    beam: PolarizedBeam,
    event: DielectricEvent,
    facet_events: u32,
    transmissions: u32,
    log_path_density: f64,
    delta: bool,
}
impl SmithDielectricSample {
    pub const fn beam(self) -> PolarizedBeam {
        self.beam
    }
    /// Net macroscopic event, determined by the final half-space.
    pub const fn event(self) -> DielectricEvent {
        self.event
    }
    pub const fn facet_events(self) -> u32 {
        self.facet_events
    }
    /// Number of microscopic transmissions, including returns to the first side.
    pub const fn transmissions(self) -> u32 {
        self.transmissions
    }
    /// A smooth/null event, or an equal-index coated path with no facet
    /// reflection. The latter is the straight-through atom of a mixed rough
    /// reflection/delta-transmission law; its hidden-path density is still not
    /// a marginal proposal probability.
    pub const fn is_delta(self) -> bool {
        self.delta
    }
    /// Joint hidden-path density in each sampled side's depth and facet-normal
    /// solid angle, with Fresnel branch/escape atoms. It is NOT an outgoing
    /// marginal direction PDF and must never be used as one in MIS.
    pub const fn log_path_density(self) -> f64 {
        self.log_path_density
    }
}

impl PolarizedBeam {
    /// Trace the complete two-sided GGX walk, using four independent variates:
    /// [free flight in (0,1), facet u in [0,1), facet v in [0,1), Fresnel branch
    /// in [0,1)]. A smooth surface consumes one quartet and uses only branch.
    /// A bare equal-index boundary consumes none and preserves the exact frame.
    ///
    /// The normal points from the physical object's inside to outside, so the
    /// original coating order reverses on each microscopic crossing. Anisotropic
    /// axes must be in the same frame; native scenes bind body axes once per hit.
    /// Every facet retains the native power Mueller weight/proposal. No G factor,
    /// albedo compensation or extra radiance eta factor multiplies this power.
    /// The caller owns medium topology after final escape. Budget exhaustion
    /// refuses the complete sample; cancellation is checked before each draw.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_smith_dielectric(
        self,
        cx: &Cx<'_>,
        normal: Vec3,
        target: impl Into<BeamMedium>,
        surface: DielectricSurface,
        settings: SmithWalkSettings,
        samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<SmithDielectricSample, SmithDielectricError> {
        self.sample_smith_dielectric_with_stack(
            cx,
            normal,
            target.into(),
            surface,
            None,
            settings,
            samples,
        )
    }

    /// The same walk with an original arbitrary spectral coating stack. A compact
    /// surface coating and separate stack cannot both own the same optical loss.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_smith_stack_dielectric(
        self,
        cx: &Cx<'_>,
        normal: Vec3,
        target: impl Into<BeamMedium>,
        surface: DielectricSurface,
        stack: DielectricCoatingStack<'_>,
        settings: SmithWalkSettings,
        samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<SmithDielectricSample, SmithDielectricError> {
        if surface.coating().is_some() {
            return Err(DielectricError::InvalidCoating.into());
        }
        self.sample_smith_dielectric_with_stack(
            cx,
            normal,
            target.into(),
            surface,
            Some(stack),
            settings,
            samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_smith_dielectric_with_stack(
        self,
        cx: &Cx<'_>,
        normal: Vec3,
        target: BeamMedium,
        surface: DielectricSurface,
        stack: Option<DielectricCoatingStack<'_>>,
        settings: SmithWalkSettings,
        mut samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<SmithDielectricSample, SmithDielectricError> {
        cx.checkpoint()?;
        let normal = unit(normal)?;
        let reverse = normal.dot(self.frame().direction()) > 0.0;
        let normal = if reverse { normal.scale(-1.0) } else { normal };
        if normal.dot(self.frame().direction()) >= 0.0 {
            return Err(BeamError::InvalidFrame.into());
        }
        let media = [self.medium(), target];
        let indices = [
            media[0].phase_index_at(self.ray().origin, self.wavelength_nm())?,
            media[1].phase_index_at(self.ray().origin, self.wavelength_nm())?,
        ];
        let interface = |reverse| match stack {
            Some(stack) => stack.interface(self.wavelength_nm(), reverse),
            None => surface.interface(self.wavelength_nm(), reverse),
        };
        let interfaces: [DielectricInterface<'_>; 2] = [interface(reverse)?, interface(!reverse)?];
        if indices[0] == indices[1] && stack.is_none() && surface.coating().is_none() {
            cx.checkpoint()?;
            return Ok(SmithDielectricSample {
                beam: Self {
                    medium: target,
                    ..self
                },
                event: DielectricEvent::Transmission,
                facet_events: 0,
                transmissions: 0,
                log_path_density: 0.0,
                delta: true,
            });
        }
        let Some(alpha) = surface.roughness_alpha() else {
            let draw = samples(0);
            validate_draw(draw)?;
            let sample = interfaces[0].sample_smooth_polarized(
                normal,
                self.frame().direction().scale(-1.0),
                indices[0],
                indices[1],
                draw[3],
            )?;
            let result =
                self.scatter_smooth_sample(normal, target, indices[0], indices[1], sample)?;
            cx.checkpoint()?;
            return Ok(SmithDielectricSample {
                beam: result.beam,
                event: result.event,
                facet_events: 1,
                transmissions: u32::from(result.event == DielectricEvent::Transmission),
                log_path_density: det::ln(result.pdf),
                delta: true,
            });
        };
        let ([alpha_x, alpha_y], tangent, bitangent) = if let Some(law) = surface.anisotropy() {
            let (tangent, bitangent) = law.tangent_frame(normal)?;
            (law.alphas(), tangent, bitangent)
        } else {
            let reference = if normal.x.abs() < 0.9 {
                Vec3::new(1.0, 0.0, 0.0)
            } else {
                Vec3::new(0.0, 1.0, 0.0)
            };
            let frame = PolarizationFrame::try_new(normal, reference)?;
            ([alpha, alpha], frame.horizontal(), frame.vertical())
        };
        let model = SmithGgx::try_new(alpha_x, alpha_y)?;
        let mut beam = self;
        let mut side = 0;
        let mut depth = 0.0;
        let mut transmissions = 0;
        let mut log_path_density = 0.0;
        for event in 0..=settings.max_facet_events() {
            cx.checkpoint()?;
            let draw = samples(event);
            validate_draw(draw)?;
            let sign = if side == 0 { 1.0 } else { -1.0 };
            let direction = beam.frame().direction();
            let local = Vec3::new(
                direction.dot(tangent),
                direction.dot(bitangent),
                direction.dot(normal),
            )
            .scale(sign);
            let step = model.sample(local, depth, [draw[0], draw[1], draw[2]])?;
            match step.event {
                SmithEvent::Escape => {
                    log_path_density += step.log_density;
                    if !log_path_density.is_finite() {
                        return Err(SmithError::NumericalRange.into());
                    }
                    cx.checkpoint()?;
                    return Ok(SmithDielectricSample {
                        beam,
                        event: if side == 0 {
                            DielectricEvent::Reflection
                        } else {
                            DielectricEvent::Transmission
                        },
                        facet_events: event,
                        transmissions,
                        log_path_density,
                        delta: indices[0] == indices[1] && transmissions == event,
                    });
                }
                SmithEvent::Reflection {
                    depth: next_depth,
                    facet,
                    ..
                } => {
                    if event == settings.max_facet_events() {
                        return Err(SmithDielectricError::FacetLimit {
                            facet_events: event,
                        });
                    }
                    let world_facet = Vec3::new(
                        tangent.x * facet.x + bitangent.x * facet.y + normal.x * facet.z,
                        tangent.y * facet.x + bitangent.y * facet.y + normal.y * facet.z,
                        tangent.z * facet.x + bitangent.z * facet.y + normal.z * facet.z,
                    )
                    .scale(sign);
                    let sample = interfaces[side].sample_smooth_polarized(
                        world_facet,
                        direction.scale(-1.0),
                        indices[side],
                        indices[1 - side],
                        draw[3],
                    )?;
                    let result = beam.scatter_smooth_sample(
                        world_facet,
                        media[1 - side],
                        indices[side],
                        indices[1 - side],
                        sample,
                    )?;
                    // SmithStep's reflected-direction Jacobian is undone:
                    // our hidden-path measure retains the sampled facet normal,
                    // then a discrete native Fresnel branch (including TIR).
                    log_path_density +=
                        step.log_density + det::ln(-4.0 * local.dot(facet)) + det::ln(result.pdf);
                    if !log_path_density.is_finite() {
                        return Err(SmithError::NumericalRange.into());
                    }
                    depth = next_depth;
                    if result.event == DielectricEvent::Transmission {
                        depth = opposite_depth(depth)?;
                        side = 1 - side;
                        transmissions += 1;
                    }
                    beam = result.beam;
                }
            }
        }
        Err(SmithDielectricError::FacetLimit {
            facet_events: settings.max_facet_events(),
        })
    }
}

fn validate_draw(draw: [f64; 4]) -> Result<(), SmithDielectricError> {
    if draw[0] <= 0.0
        || draw
            .iter()
            .any(|u| !u.is_finite() || !(0.0..1.0).contains(u))
    {
        return Err(SmithError::InvalidSample.into());
    }
    Ok(())
}

// Stable involution d -> -log(1-exp(-d)). For deep points the compensated
// log1p identity keeps the other depth when 1-exp(-d) rounds to one. Refuse
// an unrepresentable positive depth instead of silently resetting a boundary.
fn opposite_depth(depth: f64) -> Result<f64, SmithDielectricError> {
    if !depth.is_finite() || depth <= 0.0 {
        return Err(SmithError::InvalidDepth.into());
    }
    let result = if depth <= core::f64::consts::LN_2 {
        -det::ln(-det::expm1(-depth))
    } else {
        let probability = det::exp(-depth);
        let complement = 1.0 - probability;
        -det::ln(complement) + (probability - (1.0 - complement)) / complement
    };
    if !result.is_finite() || result <= 0.0 {
        return Err(SmithError::NumericalRange.into());
    }
    Ok(result)
}

#[cfg(test)]
mod tests;

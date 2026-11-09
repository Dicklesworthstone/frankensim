//! Multiple internal reflections on native rough opaque conductor beams.
//!
//! The Smith heightfield is statistical geometry, so all microscopic events
//! occur at the original geometric hit. It supplies neither metric travel,
//! sub-wavelength diffraction nor a spatial heat deposition model. The native
//! complex Fresnel/Mueller owner retains wavelength, coatings, incident medium,
//! polarization frames and phase at EVERY facet. Coherent films are local to a
//! facet; this does not introduce coherence between distinct random paths.

use fs_exec::Cx;
use fs_geom::Vec3;

use super::super::mirrors::ConductorBeamError;
use super::{BeamError, PolarizationFrame, PolarizedBeam};
use crate::conductor::microsurface::{SmithError, SmithEvent, SmithGgx, unit};
use crate::conductor::{ConductorOptics, ConductorSurface};

pub mod connection;
/// Two-sided rough glass with shared statistical heights and native refraction.
pub mod dielectric;

pub const MAX_SMITH_FACET_EVENTS: u32 = 65_536;

/// Explicit per-geometric-hit work limit. Exhaustion is an error, not a dark
/// sample, Russian roulette, or an approximation to the remaining path.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SmithWalkSettings {
    max_facet_events: u32,
}
impl SmithWalkSettings {
    pub fn try_new(max_facet_events: u32) -> Result<Self, SmithBeamError> {
        if max_facet_events == 0 || max_facet_events > MAX_SMITH_FACET_EVENTS {
            return Err(SmithBeamError::InvalidBudget);
        }
        Ok(Self { max_facet_events })
    }
    pub const fn max_facet_events(self) -> u32 {
        self.max_facet_events
    }
}

/// Failures preserve the supplied beam. Sampling callbacks are not rolled back.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SmithBeamError {
    Geometry(SmithError),
    Conductor(ConductorBeamError),
    InvalidBudget,
    FacetLimit { reflected_facets: u32 },
    Cancelled,
}
impl From<SmithError> for SmithBeamError {
    fn from(error: SmithError) -> Self {
        Self::Geometry(error)
    }
}
impl From<ConductorBeamError> for SmithBeamError {
    fn from(error: ConductorBeamError) -> Self {
        Self::Conductor(error)
    }
}
impl From<BeamError> for SmithBeamError {
    fn from(error: BeamError) -> Self {
        Self::Conductor(error.into())
    }
}
impl From<fs_exec::Cancelled> for SmithBeamError {
    fn from(_: fs_exec::Cancelled) -> Self {
        Self::Cancelled
    }
}
impl core::fmt::Display for SmithBeamError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "polarized Smith conductor: {self:?}")
    }
}
impl core::error::Error for SmithBeamError {}

/// A completed forward-power sample, deliberately separate from `BeamSample`:
/// a hidden-path density cannot serve as an outgoing marginal BSDF PDF in MIS.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmithBeamSample {
    pub(crate) beam: PolarizedBeam,
    facet_events: u32,
    log_path_density: f64,
    delta: bool,
}
impl SmithBeamSample {
    pub const fn beam(self) -> PolarizedBeam {
        self.beam
    }
    pub const fn facet_events(self) -> u32 {
        self.facet_events
    }
    /// Joint density over all sampled hidden depths and outgoing solid angles,
    /// including the final escape atom. See `SmithStep::log_density`.
    pub const fn log_path_density(self) -> f64 {
        self.log_path_density
    }
    /// True only for the native ideal mirror, whose discrete path mass is one.
    pub const fn is_delta(self) -> bool {
        self.delta
    }
}

impl PolarizedBeam {
    /// Trace all internal reflections of an isotropic Smith/GGX metal until
    /// escape, preserving polarized spectral power with the actual local medium.
    ///
    /// Each facet applies the original smooth conductor's Mueller response.
    /// No additional masking, compensation, Fresnel roulette or eta factor is
    /// applied: geometry is already sampled by the exact Smith walk.
    ///
    /// The callback receives a zero-based free-flight index. Its first draw is
    /// in (0,1), its two facet coordinates in [0,1). Escape consumes a free-flight
    /// draw but no facet. At most `max_facet_events + 1` callbacks are made so
    /// escape after the last allowed reflection can still complete. Smooth
    /// mirrors consume no draws. Cancellation is checked before each callback.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_smith_conductor(
        self,
        cx: &Cx<'_>,
        normal: Vec3,
        optics: &ConductorOptics,
        surface: ConductorSurface,
        settings: SmithWalkSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<SmithBeamSample, SmithBeamError> {
        cx.checkpoint()?;
        let mut normal = unit(normal)?;
        if normal.dot(self.frame().direction()) > 0.0 {
            normal = normal.scale(-1.0);
        }
        if normal.dot(self.frame().direction()) >= 0.0 {
            return Err(BeamError::InvalidFrame.into());
        }
        if surface.is_delta() {
            let beam = self.reflect_conductor(normal, optics)?.beam;
            cx.checkpoint()?;
            return Ok(SmithBeamSample {
                beam,
                facet_events: 1,
                log_path_density: 0.0,
                delta: true,
            });
        }
        let alpha = surface.roughness_alpha();
        let model = SmithGgx::try_new(alpha, alpha)?;
        let reference = if normal.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        let basis = PolarizationFrame::try_new(normal, reference)?;
        self.sample_smith_conductor_in_frame(
            cx,
            optics,
            model,
            [basis.horizontal(), basis.vertical(), normal],
            settings,
            samples,
        )
    }

    // Both the isotropic adapter and the original body-locked brushed owner
    // admit this right-handed orthonormal frame before entering the same walk.
    // There is no frame reconstruction, allocation, or finish lookup per facet.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn sample_smith_conductor_in_frame(
        self,
        cx: &Cx<'_>,
        optics: &ConductorOptics,
        model: SmithGgx,
        [tangent, bitangent, normal]: [Vec3; 3],
        settings: SmithWalkSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<SmithBeamSample, SmithBeamError> {
        self.sample_smith_conductor_in_frame_observed(
            cx,
            optics,
            model,
            [tangent, bitangent, normal],
            settings,
            samples,
            |_, _| Ok(()),
        )
    }

    // A local source observes the actual arriving Mueller prefix and facet
    // before reflection. It changes neither geometry, draws nor continuation.
    // The generic failure retains the source owner's typed refusal, and the
    // completed sample is published only after both work streams succeed.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn sample_smith_conductor_in_frame_observed<E: From<SmithBeamError>>(
        self,
        cx: &Cx<'_>,
        optics: &ConductorOptics,
        model: SmithGgx,
        [tangent, bitangent, normal]: [Vec3; 3],
        settings: SmithWalkSettings,
        mut samples: impl FnMut(u32) -> [f64; 3],
        mut observe: impl FnMut(PolarizedBeam, Vec3) -> Result<(), E>,
    ) -> Result<SmithBeamSample, E> {
        let mut beam = self;
        let mut depth = 0.0;
        let mut log_path_density = 0.0;
        for event in 0..=settings.max_facet_events {
            cx.checkpoint().map_err(SmithBeamError::from)?;
            let direction = beam.frame().direction();
            let local = Vec3::new(
                direction.dot(tangent),
                direction.dot(bitangent),
                direction.dot(normal),
            );
            let next = model
                .sample(local, depth, samples(event))
                .map_err(SmithBeamError::from)?;
            log_path_density += next.log_density;
            if !log_path_density.is_finite() {
                return Err(SmithBeamError::from(SmithError::NumericalRange).into());
            }
            match next.event {
                SmithEvent::Escape => {
                    cx.checkpoint().map_err(SmithBeamError::from)?;
                    return Ok(SmithBeamSample {
                        beam,
                        facet_events: event,
                        log_path_density,
                        delta: false,
                    });
                }
                SmithEvent::Reflection {
                    depth: next_depth,
                    facet,
                    ..
                } => {
                    if event == settings.max_facet_events {
                        return Err(SmithBeamError::FacetLimit {
                            reflected_facets: event,
                        }
                        .into());
                    }
                    let world_facet = Vec3::new(
                        tangent.x * facet.x + bitangent.x * facet.y + normal.x * facet.z,
                        tangent.y * facet.x + bitangent.y * facet.y + normal.y * facet.z,
                        tangent.z * facet.x + bitangent.z * facet.y + normal.z * facet.z,
                    );
                    observe(beam, world_facet)?;
                    beam = beam
                        .reflect_conductor(world_facet, optics)
                        .map_err(SmithBeamError::from)?
                        .beam;
                    depth = next_depth;
                }
            }
        }
        Err(SmithBeamError::FacetLimit {
            reflected_facets: settings.max_facet_events,
        }
        .into())
    }
}

#[cfg(test)]
mod tests;

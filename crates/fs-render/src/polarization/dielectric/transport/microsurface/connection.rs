//! Stochastic polarized fixed-direction response of the complete Smith metal.
//!
//! A sampled hidden collision contributes phase density times escape to the
//! requested light direction. The native Mueller law is evaluated at that
//! connection's facet, BEFORE the independently sampled continuation facet.
//! The first-order height integral is analytic; only higher orders are random.
//! This is the next-event estimator of Heitz et al. (2016), not a fit, an albedo
//! compensation, or a marginal PDF. Do not put a random response in a nonlinear
//! MIS heuristic. A finite facet limit refuses the whole estimate.
//!
//! Frames, coatings and reciprocal source contraction retain their original
//! owners. Heights are statistical and add no geometric distance or phase.

use super::{SmithBeamError, SmithWalkSettings};
use crate::conductor::microsurface::connection::SmithConnection;
use crate::conductor::microsurface::{SmithEvent, SmithGgx, unit};
use crate::conductor::{ConductorOptics, ConductorSurface};
use crate::instances::RigidTransform;
use crate::polarization::Stokes;
use crate::polarization::dielectric::transport::connection::{
    checked_product, normalized_camera, reciprocal_contraction,
};
use crate::polarization::dielectric::transport::{BeamError, PolarizationFrame, PolarizedBeam};
use crate::surface_detail::SurfaceDetail;
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};

/// A stochastic radiance BRDF contraction, intentionally without a PDF.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmithConnectionResponse {
    factor_per_sr: f64,
    facet_events: u32,
}
impl SmithConnectionResponse {
    /// Multiply by camera-prefix intensity, physical source irradiance on the
    /// geometric receiver plane, and the original camera's local n² ratio.
    pub const fn factor_per_sr(self) -> f64 {
        self.factor_per_sr
    }
    /// Actual continuation reflections; the first analytic connection is not
    /// an extra sampled facet. Escape consumes one additional free-flight draw.
    pub const fn facet_events(self) -> u32 {
        self.facet_events
    }
    const fn zero() -> Self {
        Self {
            factor_per_sr: 0.0,
            facet_events: 0,
        }
    }
}

impl PolarizedBeam {
    /// Evaluate all reflection orders toward one fixed source-arrival frame.
    /// The source propagates toward the surface and has unit Stokes intensity.
    /// Samples follow the original `sample_smith_conductor` contract. Their
    /// stream must be independent of any randomized light-direction proposal.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_smith_conductor_connection(
        self,
        cx: &Cx<'_>,
        normal: Vec3,
        optics: &ConductorOptics,
        surface: ConductorSurface,
        source_frame: PolarizationFrame,
        normalized_source: Stokes,
        settings: SmithWalkSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<SmithConnectionResponse, SmithBeamError> {
        cx.checkpoint()?;
        let camera = normalized_camera(self, normalized_source)?;
        let mut normal = unit(normal)?;
        if normal.dot(camera.frame().direction()) > 0.0 {
            normal = normal.scale(-1.0);
        }
        if normal.dot(camera.frame().direction()) >= 0.0 {
            return Err(BeamError::InvalidFrame.into());
        }
        if surface.is_delta() || normal.dot(source_frame.direction()) >= 0.0 {
            return Ok(SmithConnectionResponse::zero());
        }
        let reference = if normal.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        let basis = PolarizationFrame::try_new(normal, reference)?;
        let alpha = surface.roughness_alpha();
        camera.smith_connection_in_frame(
            cx,
            optics,
            SmithGgx::try_new(alpha, alpha)?,
            [basis.horizontal(), basis.vertical(), normal],
            source_frame,
            normalized_source,
            settings,
            samples,
        )
    }

    /// The same estimator in the original body-locked brushed frame, retaining
    /// anisotropic slopes, placement, roughness marks and undefined-pole fallback.
    /// An authored tilted normal retains the existing empirical radiance law:
    /// convert its shading-plane cosine to the geometric receiving cosine ONCE,
    /// with no per-facet adjoint factor or new energy-conservation claim.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_instanced_smith_brushed_conductor_connection(
        self,
        cx: &Cx<'_>,
        optics: &ConductorOptics,
        detail: &SurfaceDetail,
        transform: RigidTransform,
        local_point: Point3,
        world_normal: Vec3,
        source_frame: PolarizationFrame,
        normalized_source: Stokes,
        settings: SmithWalkSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<SmithConnectionResponse, SmithBeamError> {
        cx.checkpoint()?;
        let camera = normalized_camera(self, normalized_source)?;
        let Some(frame) =
            camera.smith_instanced_brushed_frame(detail, transform, local_point, world_normal)?
        else {
            let surface = ConductorSurface::try_rough(detail.alpha_pair().0)
                .map_err(super::ConductorBeamError::from)?;
            return camera.evaluate_smith_conductor_connection(
                cx,
                world_normal,
                optics,
                surface,
                source_frame,
                normalized_source,
                settings,
                samples,
            );
        };
        let outgoing = source_frame.direction().scale(-1.0);
        if camera
            .smith_brushed_endpoint_weight(&frame, outgoing, true)?
            .is_none()
        {
            return Ok(SmithConnectionResponse::zero());
        }
        let (ax, ay) = detail.alpha_pair_for_mark(frame.mark_covered);
        let mut result = camera.smith_connection_in_frame(
            cx,
            optics,
            SmithGgx::try_new(ax, ay)?,
            [frame.tangent, frame.bitangent, frame.shading_normal],
            source_frame,
            normalized_source,
            settings,
            samples,
        )?;
        let cosine_ratio =
            frame.shading_normal.dot(outgoing) / frame.geometric_normal.dot(outgoing);
        result.factor_per_sr = checked_product(result.factor_per_sr, cosine_ratio)?;
        cx.checkpoint()?;
        Ok(result)
    }

    #[allow(clippy::too_many_arguments)]
    fn smith_connection_in_frame(
        self,
        cx: &Cx<'_>,
        optics: &ConductorOptics,
        model: SmithGgx,
        basis: [Vec3; 3],
        source_frame: PolarizationFrame,
        normalized_source: Stokes,
        settings: SmithWalkSettings,
        mut samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<SmithConnectionResponse, SmithBeamError> {
        let outgoing = local(basis, source_frame.direction().scale(-1.0));
        let mut sum = Sum::default();
        if let Some(connection) =
            model.single_scattering_connection(local(basis, self.frame().direction()), outgoing)?
        {
            sum.add(contribution(
                self,
                optics,
                basis,
                connection,
                source_frame,
                normalized_source,
            )?)?;
        }
        let mut beam = self;
        let mut depth = 0.0;
        for event in 0..=settings.max_facet_events() {
            cx.checkpoint()?;
            let incoming = local(basis, beam.frame().direction());
            let step = model.sample(incoming, depth, samples(event))?;
            match step.event {
                SmithEvent::Escape => {
                    cx.checkpoint()?;
                    return Ok(SmithConnectionResponse {
                        factor_per_sr: sum.value()?,
                        facet_events: event,
                    });
                }
                SmithEvent::Reflection {
                    depth: next_depth,
                    facet,
                    ..
                } => {
                    if event == settings.max_facet_events() {
                        return Err(SmithBeamError::FacetLimit {
                            reflected_facets: event,
                        });
                    }
                    // First-order height is already integrated analytically.
                    // Never count its random connection a second time.
                    if event != 0 {
                        if let Some(connection) =
                            model.connection_at_depth(incoming, next_depth, outgoing)?
                        {
                            sum.add(contribution(
                                beam,
                                optics,
                                basis,
                                connection,
                                source_frame,
                                normalized_source,
                            )?)?;
                        }
                    }
                    beam = beam.reflect_conductor(world(basis, facet), optics)?.beam;
                    depth = next_depth;
                    if beam.power().i == 0.0 {
                        // Exact optical extinction cannot revive at later facets.
                        cx.checkpoint()?;
                        return Ok(SmithConnectionResponse {
                            factor_per_sr: sum.value()?,
                            facet_events: event + 1,
                        });
                    }
                }
            }
        }
        Err(SmithBeamError::FacetLimit {
            reflected_facets: settings.max_facet_events(),
        })
    }
}

fn contribution(
    beam: PolarizedBeam,
    optics: &ConductorOptics,
    basis: [Vec3; 3],
    connection: SmithConnection,
    source_frame: PolarizationFrame,
    source: Stokes,
) -> Result<f64, SmithBeamError> {
    let reflected = beam
        .reflect_conductor(world(basis, connection.facet), optics)?
        .beam;
    let optical = reciprocal_contraction(reflected, source_frame, source)?;
    Ok(checked_product(connection.factor_per_sr, optical)?)
}

fn local(basis: [Vec3; 3], direction: Vec3) -> Vec3 {
    Vec3::new(
        direction.dot(basis[0]),
        direction.dot(basis[1]),
        direction.dot(basis[2]),
    )
}
fn world(basis: [Vec3; 3], direction: Vec3) -> Vec3 {
    Vec3::new(
        basis[0].x * direction.x + basis[1].x * direction.y + basis[2].x * direction.z,
        basis[0].y * direction.x + basis[1].y * direction.y + basis[2].y * direction.z,
        basis[0].z * direction.x + basis[1].z * direction.y + basis[2].z * direction.z,
    )
}
#[derive(Default)]
struct Sum {
    total: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), SmithBeamError> {
        let next = self.total + value;
        self.correction += if self.total >= value {
            (self.total - next) + value
        } else {
            (value - next) + self.total
        };
        self.total = next;
        self.value().map(|_| ())
    }
    fn value(&self) -> Result<f64, SmithBeamError> {
        let value = self.total + self.correction;
        if value.is_finite() && value >= 0.0 {
            Ok(value)
        } else {
            Err(BeamError::NumericalRange.into())
        }
    }
}

#[cfg(test)]
mod tests;

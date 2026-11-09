//! Fixed-direction polarized response of the complete two-sided Smith glass.
//!
//! Integrate the first height analytically; at every later sampled height,
//! connect to the specified source with the native
//! facet Mueller operator, the visible-normal/refraction Jacobian, and escape
//! from the correct side of the SAME height. Continue independently with the
//! original full VNDF and native Fresnel branch. This is Heitz et al. (2016),
//! sections 6--7, including all reflection/transmission orders.
//! <https://doi.org/10.1145/2897824.2925943>
//!
//! The transmission solid-angle Jacobian is
//! n_t^2 |wi.m| / (n_i wo.m + n_t wi.m)^2; see Walter et al. (2007),
//! equation 17, <https://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.pdf>.
//! The result is a radiance BSDF, not an outgoing marginal PDF. Its net
//! (initial_index / source_index)^2 conversion is applied once, after all
//! microscopic power transport. The original camera prefix index conversion
//! remains with its radiometric owner. No stochastic value belongs in a
//! nonlinear MIS weight.
//!
//! Equal-index coated rough surfaces mix continuous scattering and a null
//! transmission atom. Their continuous component connects at the final facet
//! reflection, then samples all subsequent null transmissions. The atom remains
//! exclusively with the original path sampler, never in this per-steradian
//! response. The continuous-only capability query deliberately excludes mixed
//! laws so a caller cannot suppress their unsampled delta endpoint.

use super::{
    BeamError, BeamMedium, DielectricCoatingStack, DielectricError, DielectricEvent,
    DielectricInterface, DielectricSurface, PolarizationFrame, PolarizedBeam, SmithDielectricError,
    SmithError, SmithEvent, SmithGgx, SmithWalkSettings, opposite_depth, unit, validate_draw,
};
use crate::polarization::Stokes;
use crate::polarization::dielectric::transport::connection::{
    normalized_camera, reciprocal_contraction,
};
use crate::polarization::dielectric::transport::microsurface::connection::{PositiveTerm, Sum};
use fs_exec::Cx;
use fs_geom::Vec3;
use fs_math::det;

mod mixed;

/// The continuous all-order radiance BSDF contraction, intentionally without
/// PDF. An equal-index coated interface additionally has a straight-through
/// atom, carried only by `sample_smith_dielectric`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmithDielectricConnectionResponse {
    factor_per_sr: f64,
    event: DielectricEvent,
    facet_events: u32,
    transmissions: u32,
}
impl SmithDielectricConnectionResponse {
    /// Multiply by prefix intensity, source irradiance on the actual geometric
    /// receiver plane, and the original camera-prefix index ratio squared.
    /// The local initial/source index ratio is already included here.
    pub const fn factor_per_sr(self) -> f64 {
        self.factor_per_sr
    }
    /// Net connection event, from the original camera half-space to the source.
    pub const fn event(self) -> DielectricEvent {
        self.event
    }
    /// Actual continuation collisions, excluding any mixed-law null tails.
    pub const fn facet_events(self) -> u32 {
        self.facet_events
    }
    pub const fn transmissions(self) -> u32 {
        self.transmissions
    }
    fn zero(event: DielectricEvent) -> Self {
        Self {
            factor_per_sr: 0.0,
            event,
            facet_events: 0,
            transmissions: 0,
        }
    }
}

impl PolarizedBeam {
    /// Whether this interface's entire law is continuous. Equal-index rough
    /// films return false even though fixed-direction evaluation supplies their
    /// continuous component: a lighting partition must separately retain their
    /// straight-through atom. This query performs no random draws or optical
    /// coating solve.
    pub fn smith_dielectric_connection_is_continuous(
        self,
        opposite_medium: impl Into<BeamMedium>,
        surface: DielectricSurface,
    ) -> Result<bool, BeamError> {
        if surface.is_delta() {
            return Ok(false);
        }
        let initial = self
            .medium()
            .phase_index_at(self.ray().origin, self.wavelength_nm())?;
        let target = opposite_medium
            .into()
            .phase_index_at(self.ray().origin, self.wavelength_nm())?;
        Ok(initial != target)
    }

    /// Evaluate full rough-glass reflection or transmission toward one source.
    /// The source frame propagates FROM the source TO the receiver, with Stokes
    /// intensity exactly one. The target is the original OTHER boundary medium
    /// even for reflection. Normals and body-bound anisotropy use the native
    /// sampler's conventions. Four independent variates follow that sampler's
    /// contract and must be independent of any sampled light direction.
    ///
    /// At matched indices the response contains every path with at least one
    /// facet reflection, including later null transmissions, but excludes the
    /// straight-through atom. The quartet callback then receives a consecutive
    /// draw ordinal across both continuation and null-tail walks. The facet
    /// budget bounds their TOTAL sampled collisions, not just one path; escape
    /// draws do not consume facets. Exhaustion refuses the entire estimate.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_smith_dielectric_connection(
        self,
        cx: &Cx<'_>,
        normal: Vec3,
        opposite_medium: impl Into<BeamMedium>,
        surface: DielectricSurface,
        source_frame: PolarizationFrame,
        normalized_source: Stokes,
        settings: SmithWalkSettings,
        samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<SmithDielectricConnectionResponse, SmithDielectricError> {
        self.evaluate_smith_dielectric_connection_with_stack(
            cx,
            normal,
            opposite_medium.into(),
            surface,
            None,
            source_frame,
            normalized_source,
            settings,
            samples,
        )
    }

    /// The same estimator with an original arbitrary spectral stack, retaining
    /// its actual layer order on every microscopic crossing. A compact surface
    /// coating and a separate stack cannot both own the same optical loss.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_smith_stack_dielectric_connection(
        self,
        cx: &Cx<'_>,
        normal: Vec3,
        opposite_medium: impl Into<BeamMedium>,
        surface: DielectricSurface,
        stack: DielectricCoatingStack<'_>,
        source_frame: PolarizationFrame,
        normalized_source: Stokes,
        settings: SmithWalkSettings,
        samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<SmithDielectricConnectionResponse, SmithDielectricError> {
        if surface.coating().is_some() {
            return Err(DielectricError::InvalidCoating.into());
        }
        self.evaluate_smith_dielectric_connection_with_stack(
            cx,
            normal,
            opposite_medium.into(),
            surface,
            Some(stack),
            source_frame,
            normalized_source,
            settings,
            samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate_smith_dielectric_connection_with_stack(
        self,
        cx: &Cx<'_>,
        normal: Vec3,
        target: BeamMedium,
        surface: DielectricSurface,
        stack: Option<DielectricCoatingStack<'_>>,
        source_frame: PolarizationFrame,
        normalized_source: Stokes,
        settings: SmithWalkSettings,
        mut samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<SmithDielectricConnectionResponse, SmithDielectricError> {
        cx.checkpoint()?;
        let camera = normalized_camera(self, normalized_source)?;
        let normal = unit(normal)?;
        let reverse = normal.dot(camera.frame().direction()) > 0.0;
        let normal = if reverse { normal.scale(-1.0) } else { normal };
        if normal.dot(camera.frame().direction()) >= 0.0 {
            return Err(BeamError::InvalidFrame.into());
        }
        let source_direction = source_frame.direction().scale(-1.0);
        let exit_cosine = normal.dot(source_direction);
        let source_side = usize::from(exit_cosine < 0.0);
        let event = if source_side == 0 {
            DielectricEvent::Reflection
        } else {
            DielectricEvent::Transmission
        };
        let Some(alpha) = surface.roughness_alpha() else {
            return Ok(SmithDielectricConnectionResponse::zero(event));
        };
        let media = [camera.medium(), target];
        let indices = [
            media[0].phase_index_at(camera.ray().origin, camera.wavelength_nm())?,
            media[1].phase_index_at(camera.ray().origin, camera.wavelength_nm())?,
        ];
        let interface = |reverse| match stack {
            Some(stack) => stack.interface(camera.wavelength_nm(), reverse),
            None => surface.interface(camera.wavelength_nm(), reverse),
        };
        let interfaces: [DielectricInterface<'_>; 2] = [interface(reverse)?, interface(!reverse)?];
        if indices[0] == indices[1]
            && stack.is_none_or(|stack| stack.layer_count() == 0)
            && surface.coating().is_none()
        {
            return Ok(SmithDielectricConnectionResponse::zero(event));
        }
        if exit_cosine == 0.0 {
            return Ok(SmithDielectricConnectionResponse::zero(event));
        }
        let (alphas, tangent, bitangent) = if let Some(law) = surface.anisotropy() {
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
        let basis = [tangent, bitangent, normal];
        let model = SmithGgx::try_new(alphas[0], alphas[1])?;
        let outgoing = local(basis, source_direction);
        if indices[0] == indices[1] {
            return mixed::evaluate(
                cx,
                camera,
                media,
                indices[0],
                interfaces,
                model,
                basis,
                source_frame,
                normalized_source,
                settings,
                samples,
            );
        }
        let log_eta = 2.0 * (det::ln(indices[0]) - det::ln(indices[source_side]));
        let mut sum = Sum::default();
        let mut beam = camera;
        let mut depth = 0.0;
        let mut side = 0;
        let mut transmissions = 0;
        for collision in 0..=settings.max_facet_events() {
            cx.checkpoint()?;
            let draw = samples(collision);
            validate_draw(draw)?;
            let sign = if side == 0 { 1.0 } else { -1.0 };
            let direction = beam.frame().direction();
            let incoming = local(basis, direction).scale(sign);
            let step = model.sample(incoming, depth, [draw[0], draw[1], draw[2]])?;
            match step.event {
                SmithEvent::Escape => {
                    cx.checkpoint()?;
                    return Ok(SmithDielectricConnectionResponse {
                        factor_per_sr: sum.value()?,
                        event,
                        facet_events: collision,
                        transmissions,
                    });
                }
                SmithEvent::Reflection {
                    depth: next_depth,
                    facet,
                    ..
                } => {
                    if collision == settings.max_facet_events() {
                        return Err(SmithDielectricError::FacetLimit {
                            facet_events: collision,
                        });
                    }
                    if let Some(connection) = connection_at_depth(
                        model,
                        incoming,
                        if collision == 0 {
                            None
                        } else {
                            Some(next_depth)
                        },
                        outgoing.scale(sign),
                        indices[side],
                        indices[1 - side],
                        side == source_side,
                    )? {
                        let world_facet = world(basis, connection.facet).scale(sign);
                        let optics = interfaces[side].fresnel_polarized(
                            (-incoming.dot(connection.facet)).clamp(0.0, 1.0),
                            indices[side],
                            indices[1 - side],
                        )?;
                        let matrix = if side == source_side {
                            optics.reflection
                        } else {
                            optics.transmission
                        };
                        let connection_event = if side == source_side {
                            DielectricEvent::Reflection
                        } else {
                            DielectricEvent::Transmission
                        };
                        let connected = beam.scatter(
                            world_facet,
                            source_direction,
                            matrix,
                            connection_event,
                            media[1 - side],
                        )?;
                        let optical =
                            reciprocal_contraction(connected, source_frame, normalized_source)?;
                        sum.add(PositiveTerm::from_log_factor(
                            connection.log_factor_per_sr + log_eta,
                            optical,
                        )?)?;
                    }
                    // The connection uses the forced facet above. Only NOW
                    // apply the independently drawn continuation facet/branch.
                    let world_facet = world(basis, facet).scale(sign);
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
                    depth = next_depth;
                    if result.event == DielectricEvent::Transmission {
                        depth = opposite_depth(depth)?;
                        side = 1 - side;
                        transmissions += 1;
                    }
                    beam = result.beam;
                    if beam.power().i == 0.0 {
                        cx.checkpoint()?;
                        return Ok(SmithDielectricConnectionResponse {
                            factor_per_sr: sum.value()?,
                            event,
                            facet_events: collision + 1,
                            transmissions,
                        });
                    }
                }
            }
        }
        Err(SmithDielectricError::FacetLimit {
            facet_events: settings.max_facet_events(),
        })
    }
}

struct Connection {
    facet: Vec3,
    log_factor_per_sr: f64,
}

/// VNDF phase density times the escape probability, divided by source cosine.
/// With no depth, integrate the first collision exactly: same-side escape
/// gives a reciprocal sum and opposite-side escape gives the Beta height law.
#[allow(clippy::too_many_arguments)]
fn connection_at_depth(
    model: SmithGgx,
    incoming: Vec3,
    depth: Option<f64>,
    outgoing: Vec3,
    ni: f64,
    nt: f64,
    reflection: bool,
) -> Result<Option<Connection>, SmithDielectricError> {
    let wo = incoming.scale(-1.0);
    let co = outgoing.z.abs();
    if co == 0.0 || (reflection && outgoing.z < 0.0) || (!reflection && outgoing.z > 0.0) {
        return Ok(None);
    }
    let half = if reflection {
        Vec3::new(wo.x + outgoing.x, wo.y + outgoing.y, wo.z + outgoing.z)
    } else {
        // Stable near index matching and antiparallel directions. Avoid
        // subtracting two large almost equal index-scaled components.
        let difference = nt - ni;
        Vec3::new(
            ni * (wo.x + outgoing.x) + difference * outgoing.x,
            ni * (wo.y + outgoing.y) + difference * outgoing.y,
            ni * (wo.z + outgoing.z) + difference * outgoing.z,
        )
    };
    if half.x == 0.0 && half.y == 0.0 && half.z == 0.0 {
        return if reflection {
            Ok(None)
        } else {
            Err(SmithError::NumericalRange.into())
        };
    }
    let mut facet = unit(half)?;
    if facet.z < 0.0 {
        facet = facet.scale(-1.0);
    }
    let ci = wo.dot(facet);
    let ct = outgoing.dot(facet);
    if facet.z <= 0.0 || ci <= 0.0 || (reflection && ct <= 0.0) || (!reflection && ct >= 0.0) {
        return Ok(None);
    }
    let area = model.projected_area(wo);
    if area <= 0.0 || !area.is_finite() {
        return Err(SmithError::NumericalRange.into());
    }
    let (log_phase, exit_direction) = if reflection {
        (
            det::ln(model.ndf(facet)) - det::ln(4.0) - det::ln(area),
            outgoing,
        )
    } else {
        // half.dot(facet) equals ni*ci+nt*ct, retaining near-index contrast
        // that would be lost by cancellation in those two O(index) terms.
        let denominator = half.dot(facet).abs();
        if denominator <= 0.0 || !denominator.is_finite() {
            return Err(SmithError::NumericalRange.into());
        }
        (
            det::ln(model.ndf(facet)) + det::ln(ci) + 2.0 * det::ln(nt) + det::ln(-ct)
                - det::ln(area)
                - 2.0 * det::ln(denominator),
            outgoing.scale(-1.0),
        )
    };
    let escape_area = model.projected_area(exit_direction.scale(-1.0));
    let lambda_exit = escape_area / co;
    let log_height = if let Some(depth) = depth {
        let exit_depth = if reflection {
            depth
        } else {
            opposite_depth(depth)?
        };
        if escape_area == 0.0 || exit_depth == 0.0 {
            0.0
        } else {
            -exit_depth * lambda_exit
        }
    } else {
        if incoming.z >= 0.0 {
            return Err(SmithError::InvalidDirection.into());
        }
        let arrival_rate = 1.0 + model.projected_area(incoming) / -incoming.z;
        if reflection {
            det::ln(arrival_rate) - det::ln(arrival_rate + lambda_exit)
        } else {
            det::ln(arrival_rate) + log_beta(arrival_rate, 1.0 + lambda_exit)?
        }
    };
    let log_factor_per_sr = log_phase - det::ln(co) + log_height;
    if !log_factor_per_sr.is_finite() {
        return Err(SmithError::NumericalRange.into());
    }
    Ok(Some(Connection {
        facet,
        log_factor_per_sr,
    }))
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

// Positive-argument log Beta. Smith's two shapes are >= 1. Recurrences move
// both to >= 8 in at most fourteen steps. The combined Stirling expression
// avoids cancellation between three large log-Gammas. The first omitted
// positive-real Stirling term bounds the series remainder (DLMF 5.11(ii));
// at x=8 it is below 2e-17. Floating arithmetic is still ordinary f64.
// https://dlmf.nist.gov/5.11 and https://dlmf.nist.gov/5.12
fn log_beta(mut a: f64, mut b: f64) -> Result<f64, SmithDielectricError> {
    if !a.is_finite() || !b.is_finite() || a < 1.0 || b < 1.0 {
        return Err(SmithError::NumericalRange.into());
    }
    if a == 1.0 {
        return Ok(-det::ln(b));
    }
    if b == 1.0 {
        return Ok(-det::ln(a));
    }
    if a > b {
        core::mem::swap(&mut a, &mut b);
    }
    let mut recurrence = 0.0;
    while a < 8.0 {
        recurrence += det::ln(a + b) - det::ln(a);
        a += 1.0;
    }
    while b < 8.0 {
        recurrence += det::ln(a + b) - det::ln(b);
        b += 1.0;
    }
    if a > b {
        core::mem::swap(&mut a, &mut b);
    }
    let ratio = a / b;
    let sum = a + b;
    let plus = 1.0 + ratio;
    // Compensated log1p: preserve the small positive logarithm even when
    // 1+ratio rounds to one, before multiplication by the large other shape.
    let log_plus = det::ln(plus) + (ratio - (plus - 1.0)) / plus;
    let result =
        recurrence + (a - 0.5) * det::ln(ratio) - (sum - 1.0) * log_plus - 0.5 * det::ln(sum)
            + 0.5 * det::ln(core::f64::consts::TAU)
            + stirling_remainder(a)
            + stirling_remainder(b)
            - stirling_remainder(sum);
    if result.is_finite() {
        Ok(result)
    } else {
        Err(SmithError::NumericalRange.into())
    }
}
fn stirling_remainder(x: f64) -> f64 {
    let inverse = 1.0 / x;
    let square = inverse * inverse;
    inverse
        * (1.0 / 12.0
            + square
                * (-1.0 / 360.0
                    + square
                        * (1.0 / 1260.0
                            + square
                                * (-1.0 / 1680.0
                                    + square
                                        * (1.0 / 1188.0
                                            + square
                                                * (-691.0 / 360360.0
                                                    + square
                                                        * (1.0 / 156.0
                                                            + square
                                                                * (-3617.0 / 122400.0
                                                                    + square
                                                                        * (43867.0
                                                                            / 244188.0)))))))))
}

#[cfg(test)]
mod tests;

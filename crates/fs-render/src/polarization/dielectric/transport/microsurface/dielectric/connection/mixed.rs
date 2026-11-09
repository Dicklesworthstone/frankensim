//! Continuous part of an index-matched coated Smith interface.
//!
//! Every non-atomic path has a unique LAST reflection. Sample its arbitrary
//! prefix with the native walk, force that reflection toward the source, and
//! integrate the remaining transmission-only suffix. At matched bulk indices
//! those transmissions preserve direction, but their Mueller matrices, side,
//! film order and shared statistical height still change. Omitting the suffix
//! loses real continuous transmitted and reflected energy.
//!
//! Geometry samples the suffix without Fresnel roulette: multiply each native
//! transmission matrix directly. Thus no transmission branch probability is
//! divided out. The zero-transmission suffix uses exact escape; its first
//! height is analytic as in the unequal-index owner. All other suffixes are
//! sampled. Each path is counted at its final reflection exactly once.

use super::*;
use crate::conductor::microsurface::SmithStep;
use crate::polarization::Mueller;

#[allow(clippy::too_many_arguments)]
pub(super) fn evaluate(
    cx: &Cx<'_>,
    camera: PolarizedBeam,
    media: [BeamMedium; 2],
    index: f64,
    interfaces: [DielectricInterface<'_>; 2],
    model: SmithGgx,
    basis: [Vec3; 3],
    source_frame: PolarizationFrame,
    source: Stokes,
    settings: SmithWalkSettings,
    mut samples: impl FnMut(u32) -> [f64; 4],
) -> Result<SmithDielectricConnectionResponse, SmithDielectricError> {
    let direction_to_source = source_frame.direction().scale(-1.0);
    let outgoing = local(basis, direction_to_source);
    let source_side = usize::from(outgoing.z < 0.0);
    let event = if source_side == 0 {
        DielectricEvent::Reflection
    } else {
        DielectricEvent::Transmission
    };
    let mut work = Work {
        draws: 0,
        facets: 0,
        limit: settings.max_facet_events(),
        samples: &mut samples,
    };
    let mut sum = Sum::default();
    let mut beam = camera;
    let mut depth = 0.0;
    let mut side = 0;
    let mut facets = 0;
    let mut transmissions = 0;
    loop {
        let sign = if side == 0 { 1.0 } else { -1.0 };
        let direction = beam.frame().direction();
        let incoming = local(basis, direction).scale(sign);
        let (step, branch) = work.sample(cx, model, incoming, depth)?;
        let SmithEvent::Reflection {
            depth: next_depth,
            facet,
            ..
        } = step.event
        else {
            cx.checkpoint()?;
            return Ok(SmithDielectricConnectionResponse {
                factor_per_sr: sum.value()?,
                event,
                facet_events: facets,
                transmissions,
            });
        };
        let outgoing = outgoing.scale(sign);
        if let Some(connection) = reflection(model, incoming, outgoing)? {
            let world_facet = world(basis, connection.facet).scale(sign);
            let optics = interfaces[side].fresnel_polarized(
                (-incoming.dot(connection.facet)).clamp(0.0, 1.0),
                index,
                index,
            )?;
            let log_power = det::ln(beam.power().i);
            if let Some((connected, log_power)) = scatter_factored(
                beam,
                world_facet,
                direction_to_source,
                optics.reflection,
                DielectricEvent::Reflection,
                media[1 - side],
                log_power,
            )? {
                // No-null-tail term. Do not sample the first height twice.
                if let Some(immediate) = connection_at_depth(
                    model,
                    incoming,
                    if facets == 0 { None } else { Some(next_depth) },
                    outgoing,
                    index,
                    index,
                    true,
                )? {
                    sum.add(PositiveTerm::from_log_factor(
                        immediate.log_factor_per_sr + log_power,
                        reciprocal_contraction(connected, source_frame, source)?,
                    )?)?;
                }
                // All suffixes with at least one null transmission. In
                // particular the outgoing ray need not face this side's sky.
                if let Some((escaped, tail_log_power)) = transmission_tail(
                    cx, connected, next_depth, side, media, index, interfaces, model, basis,
                    log_power, &mut work,
                )? {
                    sum.add(PositiveTerm::from_log_factor(
                        connection.log_factor_per_sr + tail_log_power,
                        reciprocal_contraction(escaped, source_frame, source)?,
                    )?)?;
                }
            }
        }
        // Tail work consumes independent subsequent draws. The already sampled
        // prefix facet and branch retain their original distribution.
        let world_facet = world(basis, facet).scale(sign);
        let sample = interfaces[side].sample_smooth_polarized(
            world_facet,
            direction.scale(-1.0),
            index,
            index,
            branch,
        )?;
        let result =
            beam.scatter_smooth_sample(world_facet, media[1 - side], index, index, sample)?;
        facets += 1;
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
                facet_events: facets,
                transmissions,
            });
        }
    }
}

/// Reflection phase density / macro source cosine, WITHOUT escape. Both
/// macro hemispheres are possible because a later null transmission can put
/// the fixed source direction on the escaping side.
fn reflection(
    model: SmithGgx,
    incoming: Vec3,
    outgoing: Vec3,
) -> Result<Option<Connection>, SmithDielectricError> {
    let wo = incoming.scale(-1.0);
    let half = Vec3::new(wo.x + outgoing.x, wo.y + outgoing.y, wo.z + outgoing.z);
    if half.x == 0.0 && half.y == 0.0 && half.z == 0.0 {
        return Ok(None);
    }
    let mut facet = unit(half)?;
    if facet.z < 0.0 {
        facet = facet.scale(-1.0);
    }
    if facet.z <= 0.0 || wo.dot(facet) <= 0.0 || outgoing.dot(facet) <= 0.0 {
        return Ok(None);
    }
    let log_factor_per_sr = det::ln(model.ndf(facet))
        - det::ln(4.0)
        - det::ln(model.projected_area(wo))
        - det::ln(outgoing.z.abs());
    if !log_factor_per_sr.is_finite() {
        return Err(SmithError::NumericalRange.into());
    }
    Ok(Some(Connection {
        facet,
        log_factor_per_sr,
    }))
}

/// Raw Mueller transfer with its positive intensity scale retained in log
/// space. A long lossy null tail can be representable after multiplication by
/// its directional density even when its isolated ray power would underflow.
#[allow(clippy::too_many_arguments)]
fn scatter_factored(
    beam: PolarizedBeam,
    facet: Vec3,
    direction: Vec3,
    matrix: Mueller,
    event: DielectricEvent,
    target: BeamMedium,
    log_power: f64,
) -> Result<Option<(PolarizedBeam, f64)>, SmithDielectricError> {
    let mean = matrix.0[0][0];
    if mean == 0.0 {
        return Ok(None);
    }
    if mean < 0.0 || !mean.is_finite() {
        return Err(BeamError::NumericalRange.into());
    }
    let normalized = normalized_camera(beam, Stokes::unpolarized(1.0))?;
    let matrix = Mueller(matrix.0.map(|row| row.map(|v| v / mean)));
    let result = normalized.scatter(facet, direction, matrix, event, target)?;
    if result.power().i == 0.0 {
        return Ok(None);
    }
    let log_power = log_power + det::ln(mean) + det::ln(result.power().i);
    if !log_power.is_finite() {
        return Err(BeamError::NumericalRange.into());
    }
    Ok(Some((
        normalized_camera(result, Stokes::unpolarized(1.0))?,
        log_power,
    )))
}

#[allow(clippy::too_many_arguments)]
fn transmission_tail(
    cx: &Cx<'_>,
    mut beam: PolarizedBeam,
    mut depth: f64,
    mut side: usize,
    media: [BeamMedium; 2],
    index: f64,
    interfaces: [DielectricInterface<'_>; 2],
    model: SmithGgx,
    basis: [Vec3; 3],
    mut log_power: f64,
    work: &mut Work<'_, impl FnMut(u32) -> [f64; 4]>,
) -> Result<Option<(PolarizedBeam, f64)>, SmithDielectricError> {
    let mut transmitted = false;
    loop {
        let sign = if side == 0 { 1.0 } else { -1.0 };
        let direction = beam.frame().direction();
        let incoming = local(basis, direction).scale(sign);
        let (step, _) = work.sample(cx, model, incoming, depth)?;
        let SmithEvent::Reflection {
            depth: next_depth,
            facet,
            ..
        } = step.event
        else {
            return Ok(transmitted.then_some((beam, log_power)));
        };
        let optics = interfaces[side].fresnel_polarized(
            (-incoming.dot(facet)).clamp(0.0, 1.0),
            index,
            index,
        )?;
        let Some(result) = scatter_factored(
            beam,
            world(basis, facet).scale(sign),
            direction,
            optics.transmission,
            DielectricEvent::Transmission,
            media[1 - side],
            log_power,
        )?
        else {
            return Ok(None);
        };
        (beam, log_power) = result;
        side = 1 - side;
        depth = opposite_depth(next_depth)?;
        transmitted = true;
    }
}

struct Work<'a, F> {
    draws: u32,
    facets: u32,
    limit: u32,
    samples: &'a mut F,
}
impl<F: FnMut(u32) -> [f64; 4]> Work<'_, F> {
    fn sample(
        &mut self,
        cx: &Cx<'_>,
        model: SmithGgx,
        incoming: Vec3,
        depth: f64,
    ) -> Result<(SmithStep, f64), SmithDielectricError> {
        cx.checkpoint()?;
        let draw = (self.samples)(self.draws);
        // At most one suffix-escape draw per prefix facet, and one final
        // prefix escape: draws <= 2 * max_facet_events + 1 <= 131073.
        self.draws += 1;
        validate_draw(draw)?;
        let step = model.sample(incoming, depth, [draw[0], draw[1], draw[2]])?;
        if matches!(step.event, SmithEvent::Reflection { .. }) {
            if self.facets == self.limit {
                return Err(SmithDielectricError::FacetLimit {
                    facet_events: self.facets,
                });
            }
            self.facets += 1;
        }
        Ok((step, draw[3]))
    }
}

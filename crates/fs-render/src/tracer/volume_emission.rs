//! Direct sampling of the exact LTE source line integral.
//!
//! An emitting point proposes a direction, not an inverse-square point-light
//! contribution. Marginalizing radius gives q(w) = integral p_V(x+r*w) r^2 dr.
//! The radiance integral stops at the first surface; the proposal PDF includes
//! the whole grid, even source cells behind that surface. This removes the
//! point estimator's near-field singularity without clamping physical distance.
//! The source/first-surface decomposition follows the equation of transfer:
//! <https://pbr-book.org/4ed/Light_Transport_II_Volume_Rendering/The_Equation_of_Transfer>.

use super::*;
use crate::volumes::transport::scattering::PhaseFunction;

const EMISSION_SAMPLE_DOMAIN: u32 = 0x766f_6c65;

// Angular densities can be enormous for distant or very thin source cells.
// Normalize before addition/squaring, retaining ratios at both numeric extremes.
fn angular_weight(own: f64, other: f64, power: bool) -> f64 {
    if own == 0.0 {
        return 0.0;
    }
    if other == 0.0 {
        return 1.0;
    }
    let scale = own.max(other);
    let own = own / scale;
    let other = other / scale;
    if power {
        own * own / (own * own + other * other)
    } else {
        own / (own + other)
    }
}

/// Straight source connections are admitted at ambient opaque and phase vertices.
/// Dielectric vertices retain their complete source-integral random walk.
pub(super) enum Receiver<'a> {
    Surface {
        material: &'a Material,
        frame: SurfaceFrame,
        outgoing: Vec3,
        brushed: Option<&'a crate::surface_detail::ShadingFrame>,
    },
    Phase {
        phase: PhaseFunction,
        incoming_travel: Vec3,
    },
}

impl Receiver<'_> {
    fn admits(&self, direction: Vec3) -> bool {
        match self {
            Self::Surface { frame, .. } => frame.oriented.dot(direction) > 0.0,
            Self::Phase { .. } => true,
        }
    }

    fn shadow_origin(&self, point: Point3) -> Point3 {
        match self {
            Self::Surface { frame, .. } => point.offset(frame.oriented.scale(RAY_EPS)),
            Self::Phase { .. } => point,
        }
    }

    fn pdf(&self, direction: Vec3) -> Result<f64, TracerError> {
        match self {
            Self::Surface {
                material,
                frame,
                outgoing,
                brushed,
            } => Ok(if frame.oriented.dot(direction) > 0.0 {
                bsdf_pdf(material, frame.oriented, *outgoing, direction, *brushed)
            } else {
                0.0
            }),
            Self::Phase {
                phase,
                incoming_travel,
            } => Ok(phase.eval_cos(incoming_travel.dot(direction))?),
        }
    }

    fn response(&self, direction: Vec3, wavelength: f64) -> Result<f64, TracerError> {
        match self {
            Self::Surface {
                material,
                frame,
                outgoing,
                brushed,
            } => {
                let cosine = frame.oriented.dot(direction).max(0.0);
                if cosine == 0.0 {
                    return Ok(0.0);
                }
                Ok(opaque_bsdf_eval(
                    material,
                    frame.oriented,
                    *outgoing,
                    direction,
                    wavelength,
                    None,
                    *brushed,
                )? * cosine)
            }
            Self::Phase {
                phase,
                incoming_travel,
            } => Ok(phase.eval_cos(incoming_travel.dot(direction))?),
        }
    }
}

/// The source term precedes the free-flight outcome, so its competing PDF is
/// the original angular PDF, without the surface-light escape probability.
pub(super) fn continuation_weight(
    cx: &Cx<'_>,
    strategy: DirectStrategy,
    medium: PathMedium<'_, '_>,
    previous: Option<PreviousBsdf>,
    point: Point3,
    direction: Vec3,
) -> Result<f64, TracerError> {
    let Some(previous) = previous else {
        return Ok(1.0);
    };
    if !previous.volume_emission_competitor
        || previous.pdf.is_delta()
        || strategy == DirectStrategy::BsdfOnly
    {
        return Ok(1.0);
    }
    let source_pdf = medium
        .emission_grid()
        .emission_direction_pdf(cx, point, direction)?;
    if source_pdf == 0.0 {
        return Ok(1.0);
    }
    Ok(match strategy {
        DirectStrategy::Mis | DirectStrategy::PowerMis => angular_weight(
            previous.pdf.forward_solid_angle(),
            source_pdf,
            strategy == DirectStrategy::PowerMis,
        ),
        DirectStrategy::NeeOnly => 0.0,
        DirectStrategy::BsdfOnly => 1.0,
    })
}

#[allow(clippy::too_many_arguments)]
pub(super) fn direct_contribution(
    scene: &Scene,
    cx: &Cx<'_>,
    ray_time: Option<&PathTime>,
    settings: &Settings,
    depth: u32,
    medium: PathMedium<'_, '_>,
    point: Point3,
    receiver: Receiver<'_>,
    wavelengths: &[f64; PACKET],
    rng: &PathRng,
) -> Result<[f64; PACKET], TracerError> {
    let mut result = [0.0; PACKET];
    if settings.strategy == DirectStrategy::BsdfOnly
        || !medium.emission_grid().has_emission_proposal()
    {
        return Ok(result);
    }
    // A disjoint counter domain preserves flight, surface-light, phase, lens
    // and dielectric streams, including full-path Sobol's reserved dimensions.
    let variates = philox4x32_10(
        [rng.pixel, rng.sample, depth, EMISSION_SAMPLE_DOMAIN],
        rng.key,
    )
    .map(u32_unit);
    let Some(direction) = medium
        .emission_grid()
        .sample_emission_direction(cx, point, variates)?
    else {
        return Ok(result);
    };
    if !receiver.admits(direction) {
        return Ok(result);
    }
    let competing_pdf = receiver.pdf(direction)?;
    if !competing_pdf.is_finite() || competing_pdf < 0.0 {
        return Err(TracerError::InvalidInput);
    }
    let source_pdf = medium
        .emission_grid()
        .emission_direction_pdf(cx, point, direction)?;
    if !(source_pdf > 0.0 && source_pdf.is_finite()) {
        return Err(TracerError::Volume(
            crate::volumes::transport::VolumeTransportError::Unrepresentable,
        ));
    }
    let ray = Ray {
        origin: receiver.shadow_origin(point),
        dir: direction,
    };
    let intersection = intersect(scene, cx, &ray, ray_time)?;
    // Match the path integrator's physical endpoint convention, retaining
    // thin source cells between a numerical spawn offset and the surface.
    let (physical_direction, distance) =
        intersection.map_or((direction, f64::INFINITY), |intersection| {
            let displacement = intersection.hit.point.delta_from(point);
            let distance = displacement.norm();
            (
                if distance > 0.0 {
                    displacement.scale(1.0 / distance)
                } else {
                    direction
                },
                distance,
            )
        });
    let source = medium.source_radiance(cx, point, physical_direction, distance, wavelengths)?;
    let weight = if depth + 1 == settings.max_depth {
        1.0
    } else {
        match settings.strategy {
            DirectStrategy::Mis | DirectStrategy::PowerMis => angular_weight(
                source_pdf,
                competing_pdf,
                settings.strategy == DirectStrategy::PowerMis,
            ),
            DirectStrategy::NeeOnly => 1.0,
            DirectStrategy::BsdfOnly => return Err(TracerError::InvalidInput),
        }
    };
    for (lane, &wavelength) in wavelengths.iter().enumerate() {
        result[lane] =
            receiver.response(direction, wavelength)? * source[lane] / source_pdf * weight;
        if !result[lane].is_finite() || result[lane] < 0.0 {
            return Err(TracerError::InvalidInput);
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::angular_weight;

    #[test]
    fn source_mis_retains_complementary_weights_at_extreme_angular_densities() {
        for scale in [1.0e-300, 1.0, 1.0e300] {
            for (power, expected) in [(false, 1.0 / 3.0), (true, 1.0 / 5.0)] {
                let first = angular_weight(scale, 2.0 * scale, power);
                let reverse = angular_weight(2.0 * scale, scale, power);
                assert!((first - expected).abs() < 2.0e-15);
                assert!((first + reverse - 1.0).abs() < 2.0e-15);
            }
        }
        assert_eq!(angular_weight(1.0, 0.0, true), 1.0);
        assert_eq!(angular_weight(0.0, 1.0, false), 0.0);
    }
}

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

/// Straight source connections at ambient opaque/phase vertices and the ambient
/// side of a rough dielectric. Closed glass displaces the physical grid, so a
/// source connection cannot integrate the grid on an interior-facing ray.
pub(super) enum Receiver<'a> {
    Surface {
        material: &'a Material,
        frame: SurfaceFrame,
        outgoing: Vec3,
        brushed: Option<&'a crate::surface_detail::ShadingFrame>,
    },
    Dielectric {
        frame: SurfaceFrame,
        outgoing: Vec3,
        boundary: BoundaryMedia,
        alpha: f64,
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
            Self::Dielectric {
                frame, boundary, ..
            } => {
                let cosine = frame.oriented.dot(direction);
                (cosine > 0.0 && boundary.incident.is_none())
                    || (cosine < 0.0 && boundary.transmitted.is_none())
            }
            Self::Phase { .. } => true,
        }
    }

    fn shadow_origin(&self, point: Point3, direction: Vec3) -> Point3 {
        match self {
            Self::Surface { frame, .. } => point.offset(frame.oriented.scale(RAY_EPS)),
            Self::Dielectric { frame, .. } => {
                dielectric_spawn_origin(point, frame.geometric, direction)
            }
            Self::Phase { .. } => point,
        }
    }

    /// Cache wavelength-independent work once per connection. For a phase
    /// vertex this is also its response; rough glass has native spectral PDFs.
    fn common_pdf(&self, direction: Vec3) -> Result<Option<f64>, TracerError> {
        match self {
            Self::Surface {
                material,
                frame,
                outgoing,
                brushed,
            } => Ok(Some(bsdf_pdf(
                material,
                frame.oriented,
                *outgoing,
                direction,
                *brushed,
            ))),
            Self::Phase {
                phase,
                incoming_travel,
            } => Ok(Some(phase.eval_cos(incoming_travel.dot(direction))?)),
            Self::Dielectric { .. } => Ok(None),
        }
    }

    /// BSDF times the receiver cosine, paired with its sampling density.
    /// Dispersive interfaces split before Fresnel-event selection, so each
    /// wavelength competes with its own complete rough-dielectric proposal.
    fn evaluation(
        &self,
        direction: Vec3,
        wavelength: f64,
        common_pdf: Option<f64>,
    ) -> Result<(f64, f64), TracerError> {
        match self {
            Self::Surface {
                material,
                frame,
                outgoing,
                brushed,
            } => {
                let cosine = frame.oriented.dot(direction).max(0.0);
                if cosine == 0.0 {
                    return Ok((0.0, 0.0));
                }
                Ok((
                    opaque_bsdf_eval(
                        material,
                        frame.oriented,
                        *outgoing,
                        direction,
                        wavelength,
                        None,
                        *brushed,
                    )? * cosine,
                    common_pdf.ok_or(TracerError::InvalidInput)?,
                ))
            }
            Self::Dielectric {
                frame,
                outgoing,
                boundary,
                alpha,
            } => {
                let evaluation = evaluate_rough_dielectric(
                    frame.oriented,
                    *outgoing,
                    direction,
                    medium_ior(boundary.incident, wavelength)?,
                    medium_ior(boundary.transmitted, wavelength)?,
                    *alpha,
                )?;
                Ok((
                    evaluation.value * admitted_unit_dot(frame.oriented, direction).abs().min(1.0),
                    evaluation.pdf,
                ))
            }
            Self::Phase { .. } => {
                let value = common_pdf.ok_or(TracerError::InvalidInput)?;
                Ok((value, value))
            }
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
    #[cfg(feature = "photometry")] fixtures: Option<&photometric::FixtureLighting<'_>>,
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
    let common_pdf = receiver.common_pdf(direction)?;
    if common_pdf.is_some_and(|pdf| !pdf.is_finite() || pdf < 0.0) {
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
        origin: receiver.shadow_origin(point, direction),
        dir: direction,
    };
    #[cfg(feature = "photometry")]
    let intersection = if fixtures.is_some() {
        intersect_primitives(&scene.primitives, cx, &ray, ray_time, Some(f64::MAX))?
    } else {
        intersect(scene, cx, &ray, ray_time)?
    };
    #[cfg(not(feature = "photometry"))]
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
    #[cfg(feature = "photometry")]
    let (physical_direction, distance) = if let Some(fixtures) = fixtures {
        // Match calibrated continuation range without applying the point-light
        // singularity exclusion to this finite thermal source line integral.
        fixtures.admit_source_segment(
            point,
            direction,
            intersection.map(|intersection| intersection.hit.point),
        )?
    } else {
        (physical_direction, distance)
    };
    let source = medium.source_radiance(cx, point, physical_direction, distance, wavelengths)?;
    for (lane, &wavelength) in wavelengths.iter().enumerate() {
        let (response, competing_pdf) = receiver.evaluation(direction, wavelength, common_pdf)?;
        if !competing_pdf.is_finite() || competing_pdf < 0.0 {
            return Err(TracerError::InvalidInput);
        }
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
        result[lane] = response * source[lane] / source_pdf * weight;
        if !result[lane].is_finite() || result[lane] < 0.0 {
            return Err(TracerError::InvalidInput);
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dielectric::{BeerLambertAbsorption, CauchyIor};

    fn dielectric_receiver(
        incident: Option<DielectricGlass>,
        transmitted: Option<DielectricGlass>,
    ) -> Receiver<'static> {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        Receiver::Dielectric {
            frame: SurfaceFrame {
                oriented: normal,
                geometric: normal,
                entering: incident.is_none(),
            },
            outgoing: normal,
            boundary: BoundaryMedia {
                incident,
                transmitted,
                transition: MediumTransition::Exit {
                    boundary_primitive: 0,
                },
            },
            alpha: 0.3,
        }
    }

    #[test]
    fn g0_rough_source_connections_stay_on_the_ambient_side() {
        let glass = DielectricGlass::representative_crown();
        let positive = Vec3::new(0.0, 0.0, 1.0);
        let negative = positive.scale(-1.0);
        let entry = dielectric_receiver(None, Some(glass));
        let exit = dielectric_receiver(Some(glass), None);
        let nested = dielectric_receiver(Some(glass), Some(glass));
        assert!(entry.admits(positive));
        assert!(!entry.admits(negative));
        assert!(!exit.admits(positive));
        assert!(exit.admits(negative));
        assert!(!nested.admits(positive));
        assert!(!nested.admits(negative));
        for receiver in [entry, exit, nested] {
            assert!(!receiver.admits(Vec3::new(1.0, 0.0, 0.0)));
        }
        let origin = Point3::new(0.0, 0.0, 0.0);
        assert!(
            dielectric_receiver(None, Some(glass))
                .shadow_origin(origin, positive)
                .z
                > 0.0
        );
        assert!(
            dielectric_receiver(Some(glass), None)
                .shadow_origin(origin, negative)
                .z
                < 0.0
        );
    }

    #[test]
    fn g1_rough_source_lanes_match_normal_fresnel_and_radiance_eta_factors() {
        let glass = DielectricGlass::new(
            CauchyIor::try_new(1.4, 0.012, 0.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.0).unwrap(),
            GlassProvenance::Custom,
        );
        let reflection = dielectric_receiver(None, Some(glass));
        let transmission = dielectric_receiver(Some(glass), None);
        let d = 1.0 / (core::f64::consts::PI * 0.3_f64.powi(2));
        let mut previous_pdf = None;
        for wavelength in [450.0, 550.0, 650.0, 750.0] {
            let um = wavelength * 1.0e-3;
            let index = 1.4 + 0.012 / (um * um);
            let fresnel = ((index - 1.0) / (index + 1.0)).powi(2);
            // At normal incidence G=1 and the reflection Jacobian is 1/4.
            let (response, pdf) = reflection
                .evaluation(Vec3::new(0.0, 0.0, 1.0), wavelength, None)
                .unwrap();
            let expected = fresnel * d / 4.0;
            assert!((response / expected - 1.0).abs() < 2.0e-13);
            assert!((pdf / expected - 1.0).abs() < 2.0e-13);
            if let Some(previous) = previous_pdf {
                assert_ne!(pdf, previous, "dispersion requires a native lane PDF");
            }
            previous_pdf = Some(pdf);

            // The rough transmission Jacobian retains eta_t²/eta_i² while
            // radiance transport does not; an exit therefore gains index².
            let eta = 1.0 / index;
            let expected_response = (1.0 - fresnel) * d / (1.0 - eta).powi(2);
            let expected_pdf = expected_response * eta * eta;
            let (response, pdf) = transmission
                .evaluation(Vec3::new(0.0, 0.0, -1.0), wavelength, None)
                .unwrap();
            assert!((response / expected_response - 1.0).abs() < 2.0e-13);
            assert!((pdf / expected_pdf - 1.0).abs() < 2.0e-13);
            assert!((response / pdf / (index * index) - 1.0).abs() < 2.0e-13);

            // The directly sampled source and the native BSDF continuation
            // must partition each wavelength independently.
            for power in [false, true] {
                assert!(
                    (angular_weight(0.7, pdf, power) + angular_weight(pdf, 0.7, power) - 1.0).abs()
                        < 2.0e-15
                );
            }
        }
    }

    #[test]
    fn g0_null_index_rough_boundaries_keep_the_delta_source_path() {
        let glass = DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.0).unwrap(),
            GlassProvenance::Custom,
        );
        let exit = dielectric_receiver(Some(glass), None);
        assert!(exit.admits(Vec3::new(0.0, 0.0, -1.0)));
        assert_eq!(
            exit.evaluation(Vec3::new(0.0, 0.0, -1.0), 550.0, None)
                .unwrap(),
            (0.0, 0.0)
        );
        let sample = sample_rough_dielectric(
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, 1.0),
            1.0,
            1.0,
            0.3,
            0.2,
            0.4,
            0.7,
        )
        .unwrap()
        .unwrap();
        assert!(sample.delta);
        assert_eq!(sample.radiance_weight, 1.0);
    }

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

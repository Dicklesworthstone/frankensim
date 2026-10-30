//! Homogeneous participating media for calibrated photometric transport.
//!
//! Coefficients are wavelength-independent inverse metres. Under this declared
//! gray-medium approximation, photopic integration commutes with extinction and
//! scattering, so calibrated source lux can be transported directly. Spectral
//! absorption, density fields, multiple scattering, and scattering from light
//! already reflected by surfaces require a fuller transport model.
//!
//! The single-scattering integral reuses the volume module's Beer–Lambert and
//! Henyey–Greenstein kernels. It samples the phase-weighted incident field at
//! panel midpoints and integrates the homogeneous camera-path extinction within
//! each panel analytically. Successive refinements report their discrepancy as
//! a numerical estimate, never as a certified error bound. Shadow boundaries,
//! narrow beams, or unresolved source variation can evade that estimate.
//! The scene entry point supplies declared angular-profile boundaries and
//! refines their intervals independently; a bare callback supplies no such
//! support information.

use core::{f64::consts::PI, fmt};

use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::det;

use crate::{charts::Ray, volumes};

/// Angular redistribution of light, normalized over the full sphere.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PhaseFunction {
    /// Equal probability per solid angle: `1/(4π)` sr⁻¹.
    Isotropic,
    /// Henyey–Greenstein scattering with mean scattering cosine `g`.
    /// Positive `g` means forward scattering along the light's propagation.
    HenyeyGreenstein {
        /// Asymmetry in the strictly open interval `(-1,1)`.
        g: f64,
    },
    /// Unpolarized Rayleigh scattering: `3(1+cos²θ)/(16π)` sr⁻¹.
    Rayleigh,
}

/// Validated global, homogeneous, gray participating medium.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HomogeneousMedium {
    absorption_per_m: f64,
    scattering_per_m: f64,
    extinction_per_m: f64,
    phase: PhaseFunction,
}

impl HomogeneousMedium {
    /// Admit finite nonnegative coefficients in m⁻¹ and a normalized phase law.
    /// These values are declared model inputs, not a claim of measurement.
    pub fn try_new(
        absorption_per_m: f64,
        scattering_per_m: f64,
        phase: PhaseFunction,
    ) -> Result<Self, MediumError> {
        if !absorption_per_m.is_finite() || absorption_per_m < 0.0 {
            return Err(MediumError::Absorption);
        }
        if !scattering_per_m.is_finite() || scattering_per_m < 0.0 {
            return Err(MediumError::Scattering);
        }
        let extinction_per_m = absorption_per_m + scattering_per_m;
        if !extinction_per_m.is_finite() {
            return Err(MediumError::ExtinctionOverflow);
        }
        if let PhaseFunction::HenyeyGreenstein { g } = phase
            && (!g.is_finite() || g.abs() >= 1.0)
        {
            return Err(MediumError::Asymmetry);
        }
        Ok(Self {
            absorption_per_m,
            scattering_per_m,
            extinction_per_m,
            phase,
        })
    }

    /// No extinction or scattering.
    #[must_use]
    pub const fn vacuum() -> Self {
        Self {
            absorption_per_m: 0.0,
            scattering_per_m: 0.0,
            extinction_per_m: 0.0,
            phase: PhaseFunction::Isotropic,
        }
    }

    /// Absorption coefficient, m⁻¹.
    #[must_use]
    pub const fn absorption_per_m(self) -> f64 {
        self.absorption_per_m
    }

    /// Scattering coefficient, m⁻¹.
    #[must_use]
    pub const fn scattering_per_m(self) -> f64 {
        self.scattering_per_m
    }

    /// Total extinction `sigma_a + sigma_s`, m⁻¹.
    #[must_use]
    pub const fn extinction_per_m(self) -> f64 {
        self.extinction_per_m
    }

    /// Declared scattering law.
    #[must_use]
    pub const fn phase(self) -> PhaseFunction {
        self.phase
    }

    /// Direct, unscattered transmission over a finite nonnegative metric length.
    /// The exponent accounts for scattering out of the beam as well as absorption.
    pub fn transmittance(self, distance_m: f64) -> Result<f64, MediumError> {
        if !distance_m.is_finite() || distance_m < 0.0 {
            return Err(MediumError::Distance);
        }
        Ok(volumes::beer_lambert(self.extinction_per_m, distance_m))
    }

    /// Phase density in sr⁻¹ at the cosine between incoming and outgoing
    /// **propagation** directions. For a point light, incoming propagation is
    /// from the light to the scattering point; outgoing is toward the camera.
    /// Tiny unit-vector dot-product roundoff beyond `[-1,1]` is clamped.
    pub fn phase_per_sr(self, cos_scatter: f64) -> Result<f64, MediumError> {
        if !cos_scatter.is_finite() || cos_scatter.abs() > 1.0 + 32.0 * f64::EPSILON {
            return Err(MediumError::ScatteringCosine);
        }
        let cosine = cos_scatter.clamp(-1.0, 1.0);
        Ok(match self.phase {
            PhaseFunction::Isotropic => 1.0 / (4.0 * PI),
            // volumes::hg_pdf_cos integrates to one over d(cosθ), so restore
            // the azimuthal 1/(2π) to obtain the solid-angle phase density.
            PhaseFunction::HenyeyGreenstein { g } => volumes::hg_pdf_cos(g, cosine) / (2.0 * PI),
            PhaseFunction::Rayleigh => 3.0 * (1.0 + cosine * cosine) / (16.0 * PI),
        })
    }
}

/// Invalid homogeneous-medium input.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MediumError {
    /// Absorption must be finite and nonnegative.
    Absorption,
    /// Scattering must be finite and nonnegative.
    Scattering,
    /// The sum of the admitted coefficients cannot be represented.
    ExtinctionOverflow,
    /// Henyey–Greenstein asymmetry must be finite and strictly between −1 and 1.
    Asymmetry,
    /// Distance must be finite and nonnegative, in metres.
    Distance,
    /// The scattering cosine is nonfinite or outside the unit-vector domain.
    ScatteringCosine,
}

impl fmt::Display for MediumError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for MediumError {}

/// Explicit work and accuracy settings for the single-scattering integral.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScatterSettings {
    /// Panels in the first midpoint evaluation; must be positive.
    pub initial_panels: u32,
    /// Largest panel count in a refinement, at least `initial_panels`.
    /// Scene-level angular partitions share this total, rather than each
    /// receiving a separate budget. Every known interval needs at least one.
    /// Total callback evaluations are less than four times this limit.
    pub max_panels: u32,
    /// Requested absolute discrepancy between refinements, in cd/m².
    pub absolute_tolerance_cd_m2: f64,
    /// Requested relative discrepancy between refinements.
    pub relative_tolerance: f64,
}

impl Default for ScatterSettings {
    fn default() -> Self {
        Self {
            initial_panels: 32,
            max_panels: 4096,
            absolute_tolerance_cd_m2: 1.0e-6,
            relative_tolerance: 1.0e-4,
        }
    }
}

impl ScatterSettings {
    pub(crate) fn is_valid(self) -> bool {
        self.initial_panels > 0
            && self.max_panels >= self.initial_panels
            && self.absolute_tolerance_cd_m2.is_finite()
            && self.absolute_tolerance_cd_m2 >= 0.0
            && self.relative_tolerance.is_finite()
            && self.relative_tolerance >= 0.0
    }
}

/// Photopic single-scattered radiance and its explicitly limited uncertainty.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScatterEstimate {
    /// Single-scattered photopic radiance, cd/m²; no tone map or exposure.
    pub luminance_cd_m2: f64,
    /// Absolute difference between the last two quadratures, in cd/m².
    /// For scene-level angular partitions this is the sum of each interval's
    /// discrepancy after camera transmission, with no cancellation of errors.
    /// `None` means no refinement comparison was possible. This is a heuristic
    /// numerical estimate, not an enclosure or a bound on unmodeled physics.
    pub estimated_quadrature_error_cd_m2: Option<f64>,
    /// Incident-field callback evaluations across all refinements.
    pub evaluations: u64,
    /// Panels used in the returned quadrature; zero for an analytic zero.
    pub panels: u32,
    /// Whether the **estimated** refinement discrepancy met the requested
    /// tolerance. This does not certify resolved shadows, beams, or convergence.
    pub tolerance_met: bool,
    /// No general multiple-scattering bound is available: `None` for scattering
    /// media. An analytically zero scattering path has the exact model bound 0.
    pub omitted_multiple_scattering_bound_cd_m2: Option<f64>,
}

/// Admission, cancellation, or incident-field failure; no partial result escapes.
#[derive(Debug, PartialEq)]
pub enum ScatterError<E> {
    /// Invalid metric distance or medium input.
    Medium(MediumError),
    /// Panel counts or nonnegative finite tolerances were invalid.
    Settings,
    /// The ray origin/direction/endpoint is nonfinite, or direction is not unit.
    Ray,
    /// The incident callback returned a negative or nonfinite value.
    IncidentValue,
    /// Quadrature length, weight, or accumulated luminance was unrepresentable.
    Arithmetic,
    /// The execution context requested cancellation.
    Cancelled,
    /// The scene-specific callback refused to evaluate the incident field.
    Incident(E),
}

impl<E: fmt::Debug> fmt::Display for ScatterError<E> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl<E: fmt::Debug> std::error::Error for ScatterError<E> {}

/// Integrate `∫ sigma_s * T_camera(s) * incident(s) ds` over a camera ray.
///
/// The ray starts at the camera (or at the start of the visible medium segment),
/// has a unit direction, and extends `distance_m` metres. The callback receives
/// the world point and unit direction **toward the camera**. It must return the
/// sum of phase-weighted direct incident illuminances, in lux/sr:
/// `Σ I_v/r² * T_light * visibility * phase(cos_scatter)`. There is no surface
/// cosine at a scattering event. Medium extinction on the source leg belongs
/// to the callback; extinction on the camera leg belongs to this integrator.
///
/// Midpoint incident values use exact analytic Beer–Lambert panel weights,
/// preserving a constant incident field even in optically thick media. Panel
/// counts double (clamped to the cap) and accumulation follows fixed order with
/// compensation. The callback must describe one immutable scene throughout.
/// Cancellation is checked at each callback and before result publication.
/// No RNG, window, GPU, scene mutation, or per-panel allocation is involved.
pub fn single_scattered_luminance<E>(
    medium: &HomogeneousMedium,
    ray: &Ray,
    distance_m: f64,
    settings: ScatterSettings,
    cx: &Cx<'_>,
    mut incident: impl FnMut(Point3, Vec3) -> Result<f64, E>,
) -> Result<ScatterEstimate, ScatterError<E>> {
    cx.checkpoint().map_err(|_| ScatterError::Cancelled)?;
    medium
        .transmittance(distance_m)
        .map_err(ScatterError::Medium)?;
    if !settings.is_valid() {
        return Err(ScatterError::Settings);
    }
    let endpoint = ray.at(distance_m);
    if [
        ray.origin.x,
        ray.origin.y,
        ray.origin.z,
        ray.dir.x,
        ray.dir.y,
        ray.dir.z,
        endpoint.x,
        endpoint.y,
        endpoint.z,
    ]
    .iter()
    .any(|value| !value.is_finite())
        || (ray.dir.dot(ray.dir) - 1.0).abs() > 1.0e-12
    {
        return Err(ScatterError::Ray);
    }
    if distance_m == 0.0 || medium.scattering_per_m == 0.0 {
        cx.checkpoint().map_err(|_| ScatterError::Cancelled)?;
        return Ok(ScatterEstimate {
            luminance_cd_m2: 0.0,
            estimated_quadrature_error_cd_m2: Some(0.0),
            evaluations: 0,
            panels: 0,
            tolerance_met: true,
            omitted_multiple_scattering_bound_cd_m2: Some(0.0),
        });
    }
    let mut panels = settings.initial_panels;
    let mut previous: Option<f64> = None;
    let mut evaluations = 0_u64;
    loop {
        let width_m = distance_m / f64::from(panels);
        if width_m <= 0.0 {
            return Err(ScatterError::Arithmetic);
        }
        let panel_extinction = medium.extinction_per_m * width_m;
        let panel_weight =
            (medium.scattering_per_m / medium.extinction_per_m) * -det::expm1(-panel_extinction);
        if !panel_weight.is_finite() || panel_weight <= 0.0 {
            return Err(ScatterError::Arithmetic);
        }
        let (mut sum, mut compensation) = (0.0_f64, 0.0_f64);
        for panel in 0..panels {
            cx.checkpoint().map_err(|_| ScatterError::Cancelled)?;
            let start_m = f64::from(panel) * width_m;
            let midpoint_m = (f64::from(panel) + 0.5) * width_m;
            let value = incident(ray.at(midpoint_m), ray.dir.scale(-1.0))
                .map_err(ScatterError::Incident)?;
            if !value.is_finite() || value < 0.0 {
                return Err(ScatterError::IncidentValue);
            }
            evaluations += 1;
            let contribution =
                (value * panel_weight) * volumes::beer_lambert(medium.extinction_per_m, start_m);
            let adjusted = contribution - compensation;
            let next = sum + adjusted;
            compensation = (next - sum) - adjusted;
            sum = next;
            if !sum.is_finite() {
                return Err(ScatterError::Arithmetic);
            }
        }
        let discrepancy = previous.map(|value| (sum - value).abs());
        let tolerance_met = discrepancy.is_some_and(|error| {
            error <= settings.absolute_tolerance_cd_m2
                || error / sum.max(f64::MIN_POSITIVE) <= settings.relative_tolerance
        });
        cx.checkpoint().map_err(|_| ScatterError::Cancelled)?;
        if tolerance_met || panels == settings.max_panels {
            return Ok(ScatterEstimate {
                luminance_cd_m2: sum,
                estimated_quadrature_error_cd_m2: discrepancy,
                evaluations,
                panels,
                tolerance_met,
                omitted_multiple_scattering_bound_cd_m2: None,
            });
        }
        previous = Some(sum);
        panels = panels.saturating_mul(2).min(settings.max_panels);
    }
}

/// Integrate each declared angular interval with its own actual refinements.
/// Knots are sorted, distinct normalized ray parameters strictly inside (0,1).
/// Even allocation gives narrow angular features a useful share of the work;
/// allocation by metric interval length would starve precisely those features.
/// The public callback-only integrator remains unchanged when there are no knots.
pub(crate) fn single_scattered_luminance_partitioned<E>(
    medium: &HomogeneousMedium,
    ray: &Ray,
    distance_m: f64,
    settings: ScatterSettings,
    cx: &Cx<'_>,
    knots: &[f64],
    mut incident: impl FnMut(Point3, Vec3) -> Result<f64, E>,
) -> Result<ScatterEstimate, ScatterError<E>> {
    cx.checkpoint().map_err(|_| ScatterError::Cancelled)?;
    if !settings.is_valid()
        || knots
            .iter()
            .any(|u| !u.is_finite() || *u <= 0.0 || *u >= 1.0)
        || knots.windows(2).any(|pair| pair[0] >= pair[1])
    {
        return Err(ScatterError::Settings);
    }
    if knots.is_empty() || distance_m == 0.0 || medium.scattering_per_m == 0.0 {
        return single_scattered_luminance(medium, ray, distance_m, settings, cx, incident);
    }
    medium
        .transmittance(distance_m)
        .map_err(ScatterError::Medium)?;
    let segments = u32::try_from(knots.len())
        .ok()
        .and_then(|n| n.checked_add(1))
        .filter(|n| *n <= settings.max_panels)
        .ok_or(ScatterError::Settings)?;
    let initial_total = settings.initial_panels.max(segments);
    let share = |total: u32, index: u32| total / segments + u32::from(index < total % segments);
    let mut sums = [0.0_f64; 2];
    let mut compensation = [0.0_f64; 2];
    let mut compared_every_segment = true;
    let mut evaluations = 0_u64;
    let mut panels = 0_u32;
    let mut start_m = 0.0;
    for index in 0..segments {
        cx.checkpoint().map_err(|_| ScatterError::Cancelled)?;
        let end_m = knots
            .get(index as usize)
            .map_or(distance_m, |u| distance_m * u);
        let length_m = end_m - start_m;
        if !length_m.is_finite() || length_m <= 0.0 {
            return Err(ScatterError::Arithmetic);
        }
        let segment_ray = Ray {
            origin: ray.at(start_m),
            dir: ray.dir,
        };
        let result = single_scattered_luminance(
            medium,
            &segment_ray,
            length_m,
            ScatterSettings {
                initial_panels: share(initial_total, index),
                max_panels: share(settings.max_panels, index),
                absolute_tolerance_cd_m2: settings.absolute_tolerance_cd_m2 / f64::from(segments),
                relative_tolerance: settings.relative_tolerance,
            },
            cx,
            &mut incident,
        )?;
        let transmission = medium
            .transmittance(start_m)
            .map_err(ScatterError::Medium)?;
        compared_every_segment &= result.estimated_quadrature_error_cd_m2.is_some();
        for (channel, value) in [
            result.luminance_cd_m2 * transmission,
            result.estimated_quadrature_error_cd_m2.unwrap_or(0.0) * transmission,
        ]
        .into_iter()
        .enumerate()
        {
            let adjusted = value - compensation[channel];
            let next = sums[channel] + adjusted;
            compensation[channel] = (next - sums[channel]) - adjusted;
            sums[channel] = next;
            if !next.is_finite() || next < 0.0 {
                return Err(ScatterError::Arithmetic);
            }
        }
        evaluations += result.evaluations;
        panels += result.panels;
        start_m = end_m;
    }
    let discrepancy = compared_every_segment.then_some(sums[1]);
    let tolerance_met = discrepancy.is_some_and(|error| {
        error <= settings.absolute_tolerance_cd_m2
            || error / sums[0].max(f64::MIN_POSITIVE) <= settings.relative_tolerance
    });
    cx.checkpoint().map_err(|_| ScatterError::Cancelled)?;
    Ok(ScatterEstimate {
        luminance_cd_m2: sums[0],
        estimated_quadrature_error_cd_m2: discrepancy,
        evaluations,
        panels,
        tolerance_met,
        omitted_multiple_scattering_bound_cd_m2: None,
    })
}

#[cfg(test)]
mod tests {
    use core::convert::Infallible;

    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    use super::*;

    fn with_cx<R>(operation: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 3,
                    kernel_id: 0x5048,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            operation(&gate, &cx)
        })
    }

    fn ray() -> Ray {
        Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        }
    }

    #[test]
    fn g0_medium_admission_and_beer_lambert_composition() {
        for invalid in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(HomogeneousMedium::try_new(invalid, 0.0, PhaseFunction::Isotropic).is_err());
            assert!(HomogeneousMedium::try_new(0.0, invalid, PhaseFunction::Isotropic).is_err());
        }
        assert!(HomogeneousMedium::try_new(f64::MAX, f64::MAX, PhaseFunction::Isotropic).is_err());
        for g in [-1.0, 1.0, f64::NAN, f64::INFINITY] {
            assert!(
                HomogeneousMedium::try_new(0.1, 0.2, PhaseFunction::HenyeyGreenstein { g })
                    .is_err()
            );
        }
        let medium = HomogeneousMedium::try_new(0.1, 0.2, PhaseFunction::Isotropic).unwrap();
        let t3 = medium.transmittance(3.0).unwrap();
        let t7 = medium.transmittance(7.0).unwrap();
        assert!((t3 - (-0.9_f64).exp()).abs() < 2.0e-16);
        assert!((medium.transmittance(10.0).unwrap() - t3 * t7).abs() < 2.0e-16);
        assert_eq!(HomogeneousMedium::vacuum().transmittance(1.0e200), Ok(1.0));
        assert_eq!(medium.transmittance(0.0), Ok(1.0));
        assert!(medium.transmittance(-0.1).is_err());
        assert!(medium.transmittance(f64::INFINITY).is_err());
    }

    #[test]
    fn g0_phase_normalization_forward_sign_and_narrow_peak_stability() {
        let panels = 100_000;
        for phase in [
            PhaseFunction::Isotropic,
            PhaseFunction::HenyeyGreenstein { g: 0.7 },
            PhaseFunction::HenyeyGreenstein { g: -0.7 },
            PhaseFunction::Rayleigh,
        ] {
            let medium = HomogeneousMedium::try_new(0.0, 0.1, phase).unwrap();
            let mut integral = 0.0;
            for panel in 0..panels {
                let cosine = -1.0 + 2.0 * (f64::from(panel) + 0.5) / f64::from(panels);
                integral += medium.phase_per_sr(cosine).unwrap() * 4.0 * PI / f64::from(panels);
            }
            assert!((integral - 1.0).abs() < 2.0e-8, "{phase:?}: {integral}");
        }
        let forward =
            HomogeneousMedium::try_new(0.0, 0.1, PhaseFunction::HenyeyGreenstein { g: 0.7 })
                .unwrap();
        assert!(forward.phase_per_sr(1.0).unwrap() > forward.phase_per_sr(-1.0).unwrap());
        let narrow_g = 1.0 - f64::EPSILON;
        let narrow =
            HomogeneousMedium::try_new(0.0, 0.1, PhaseFunction::HenyeyGreenstein { g: narrow_g })
                .unwrap();
        let expected = (1.0 + narrow_g) / (4.0 * PI * (1.0 - narrow_g).powi(2));
        assert!((narrow.phase_per_sr(1.0).unwrap() / expected - 1.0).abs() < 4.0e-16);
        assert!(narrow.phase_per_sr(1.01).is_err());
    }

    #[test]
    fn g1_constant_incident_field_matches_analytic_homogeneous_slab() {
        with_cx(|_, cx| {
            for distance in [1.0e-9, 1.0, 1000.0] {
                let medium =
                    HomogeneousMedium::try_new(0.05, 0.15, PhaseFunction::Isotropic).unwrap();
                let incident = 60.0;
                let result = single_scattered_luminance(
                    &medium,
                    &ray(),
                    distance,
                    ScatterSettings::default(),
                    cx,
                    |_, towards_camera| {
                        assert_eq!(towards_camera, Vec3::new(0.0, 0.0, -1.0));
                        Ok::<_, Infallible>(incident)
                    },
                )
                .unwrap();
                let expected = incident * 0.75 * -(-0.2 * distance).exp_m1();
                assert!(
                    (result.luminance_cd_m2 / expected - 1.0).abs() < 2.0e-14,
                    "L={distance}: {result:?}, expected {expected}"
                );
                assert!(result.tolerance_met);
                assert_eq!(result.omitted_multiple_scattering_bound_cd_m2, None);
            }
        });
    }

    #[test]
    fn g1_collinear_point_source_matches_closed_form_with_both_extinction_legs() {
        with_cx(|_, cx| {
            let medium =
                HomogeneousMedium::try_new(0.03, 0.02, PhaseFunction::HenyeyGreenstein { g: 0.7 })
                    .unwrap();
            let (source_distance, camera_distance, intensity) = (10.0, 4.0, 1000.0);
            let result = single_scattered_luminance(
                &medium,
                &ray(),
                camera_distance,
                ScatterSettings {
                    initial_panels: 16,
                    max_panels: 8192,
                    absolute_tolerance_cd_m2: 1.0e-10,
                    relative_tolerance: 1.0e-8,
                },
                cx,
                |point, towards_camera| {
                    let light_distance = source_distance - point.z;
                    let propagation = Vec3::new(0.0, 0.0, -1.0);
                    Ok::<_, Infallible>(
                        intensity / (light_distance * light_distance)
                            * medium.transmittance(light_distance).unwrap()
                            * medium
                                .phase_per_sr(propagation.dot(towards_camera))
                                .unwrap(),
                    )
                },
            )
            .unwrap();
            // Along this forward-scattering path, r_light + r_camera = D.
            // ∫₀ᴸ ds/(D-s)² = 1/(D-L) - 1/D.
            let expected = intensity
                * medium.scattering_per_m()
                * medium.phase_per_sr(1.0).unwrap()
                * medium.transmittance(source_distance).unwrap()
                * (1.0 / (source_distance - camera_distance) - 1.0 / source_distance);
            assert!((result.luminance_cd_m2 / expected - 1.0).abs() < 1.0e-8);
            assert!(result.tolerance_met);
            assert!(
                (result.luminance_cd_m2 - expected).abs()
                    < result.estimated_quadrature_error_cd_m2.unwrap()
            );
        });
    }

    #[test]
    fn g1_refinement_has_second_order_and_does_not_invent_a_bound() {
        with_cx(|_, cx| {
            let medium = HomogeneousMedium::try_new(0.1, 0.2, PhaseFunction::Isotropic).unwrap();
            let integrate = |panels| {
                single_scattered_luminance(
                    &medium,
                    &ray(),
                    1.0,
                    ScatterSettings {
                        initial_panels: panels,
                        max_panels: panels,
                        absolute_tolerance_cd_m2: 0.0,
                        relative_tolerance: 0.0,
                    },
                    cx,
                    |point, _| Ok::<_, Infallible>(det::exp(0.3 * point.z)),
                )
                .unwrap()
            };
            let coarse = integrate(16);
            let fine = integrate(32);
            let coarse_error = (coarse.luminance_cd_m2 - 0.2).abs();
            let fine_error = (fine.luminance_cd_m2 - 0.2).abs();
            assert!((coarse_error / fine_error - 4.0).abs() < 1.0e-3);
            assert_eq!(coarse.estimated_quadrature_error_cd_m2, None);
            assert_eq!(coarse.omitted_multiple_scattering_bound_cd_m2, None);
            assert!(!coarse.tolerance_met);
            assert_eq!(coarse.evaluations, 16);
            assert_eq!(fine, integrate(32));
        });
    }

    #[test]
    fn g1_partitioned_slab_preserves_camera_extinction_and_total_work_budget() {
        with_cx(|_, cx| {
            let medium = HomogeneousMedium::try_new(0.1, 0.2, PhaseFunction::Isotropic).unwrap();
            let settings = ScatterSettings {
                initial_panels: 2,
                max_panels: 64,
                absolute_tolerance_cd_m2: 1.0e-12,
                relative_tolerance: 1.0e-12,
            };
            let result = single_scattered_luminance_partitioned(
                &medium,
                &ray(),
                10.0,
                settings,
                cx,
                &[0.01, 0.8],
                |_, _| Ok::<_, Infallible>(40.0),
            )
            .unwrap();
            let expected = 40.0 * (0.2 / 0.3) * -(-3.0_f64).exp_m1();
            assert!((result.luminance_cd_m2 / expected - 1.0).abs() < 2.0e-14);
            assert!(result.tolerance_met);
            assert!(result.panels <= settings.max_panels);
            assert!(result.evaluations < 4 * u64::from(settings.max_panels));
            assert_eq!(result.omitted_multiple_scattering_bound_cd_m2, None);
        });
    }

    #[test]
    fn g1_declared_support_is_refined_and_one_panel_is_not_a_comparison() {
        with_cx(|_, cx| {
            let medium = HomogeneousMedium::try_new(0.05, 0.1, PhaseFunction::Isotropic).unwrap();
            let settings = ScatterSettings {
                initial_panels: 2,
                max_panels: 128,
                absolute_tolerance_cd_m2: 0.0,
                relative_tolerance: 1.0e-12,
            };
            let field = |p: Point3, _: Vec3| {
                Ok::<_, Infallible>(if (0.05..0.07).contains(&p.z) {
                    100.0
                } else {
                    0.0
                })
            };
            let result = single_scattered_luminance_partitioned(
                &medium,
                &ray(),
                1.0,
                settings,
                cx,
                &[0.05, 0.07],
                field,
            )
            .unwrap();
            let expected =
                100.0 * (0.1 / 0.15) * (-0.15_f64 * 0.05).exp() * -(-0.15_f64 * 0.02).exp_m1();
            assert!((result.luminance_cd_m2 / expected - 1.0).abs() < 2.0e-14);
            assert!(result.tolerance_met);
            let unrefined = single_scattered_luminance_partitioned(
                &medium,
                &ray(),
                1.0,
                ScatterSettings {
                    initial_panels: 1,
                    max_panels: 3,
                    ..settings
                },
                cx,
                &[0.05, 0.07],
                field,
            )
            .unwrap();
            assert_eq!(unrefined.panels, 3);
            assert_eq!(unrefined.evaluations, 3);
            assert_eq!(unrefined.estimated_quadrature_error_cd_m2, None);
            assert!(!unrefined.tolerance_met);
            for knots in [&[-0.1][..], &[1.0][..], &[0.2, 0.2][..], &[0.3, 0.2][..]] {
                assert_eq!(
                    single_scattered_luminance_partitioned(
                        &medium,
                        &ray(),
                        1.0,
                        settings,
                        cx,
                        knots,
                        field,
                    ),
                    Err(ScatterError::Settings)
                );
            }
        });
    }

    #[test]
    fn g0_analytic_zero_refusals_and_callback_failure() {
        with_cx(|_, cx| {
            let medium = HomogeneousMedium::try_new(0.1, 0.2, PhaseFunction::Isotropic).unwrap();
            let settings = ScatterSettings::default();
            for (medium, distance) in [(HomogeneousMedium::vacuum(), 3.0), (medium, 0.0)] {
                let result = single_scattered_luminance(
                    &medium,
                    &ray(),
                    distance,
                    settings,
                    cx,
                    |_, _| -> Result<f64, Infallible> { panic!("analytic zero must not sample") },
                )
                .unwrap();
                assert_eq!(result.luminance_cd_m2, 0.0);
                assert_eq!(result.omitted_multiple_scattering_bound_cd_m2, Some(0.0));
            }
            assert_eq!(
                single_scattered_luminance(&medium, &ray(), 1.0, settings, cx, |_, _| Err(
                    "occlusion"
                )),
                Err(ScatterError::Incident("occlusion"))
            );
            for value in [-1.0, f64::NAN, f64::INFINITY] {
                assert_eq!(
                    single_scattered_luminance(&medium, &ray(), 1.0, settings, cx, |_, _| {
                        Ok::<_, Infallible>(value)
                    }),
                    Err(ScatterError::IncidentValue)
                );
            }
            for invalid in [
                ScatterSettings {
                    initial_panels: 0,
                    ..settings
                },
                ScatterSettings {
                    max_panels: 1,
                    ..settings
                },
                ScatterSettings {
                    relative_tolerance: -1.0,
                    ..settings
                },
                ScatterSettings {
                    absolute_tolerance_cd_m2: f64::NAN,
                    ..settings
                },
            ] {
                assert_eq!(
                    single_scattered_luminance(&medium, &ray(), 1.0, invalid, cx, |_, _| {
                        Ok::<_, Infallible>(1.0)
                    }),
                    Err(ScatterError::Settings)
                );
            }
            let bad_ray = Ray {
                dir: Vec3::new(0.0, 0.0, 2.0),
                ..ray()
            };
            assert_eq!(
                single_scattered_luminance(&medium, &bad_ray, 1.0, settings, cx, |_, _| {
                    Ok::<_, Infallible>(1.0)
                }),
                Err(ScatterError::Ray)
            );
        });
    }

    #[test]
    fn g4_cancellation_inside_callback_prevents_partial_publication() {
        with_cx(|gate, cx| {
            let medium = HomogeneousMedium::try_new(0.1, 0.2, PhaseFunction::Isotropic).unwrap();
            let mut calls = 0;
            let result = single_scattered_luminance(
                &medium,
                &ray(),
                1.0,
                ScatterSettings::default(),
                cx,
                |_, _| {
                    calls += 1;
                    gate.request();
                    Ok::<_, Infallible>(10.0)
                },
            );
            assert_eq!(result, Err(ScatterError::Cancelled));
            assert_eq!(calls, 1);
        });
    }
}

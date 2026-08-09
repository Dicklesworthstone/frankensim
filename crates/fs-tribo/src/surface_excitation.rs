//! Contact-patch-filtered excitation from named surface-height traces.
//!
//! This module is deliberately upstream of acoustics. It converts actual
//! relative contact-path kinematics plus a measured/declared spatial texture
//! into a linearized normal-force perturbation. A structural solver may then
//! consume that force; no sound pressure, loudness target, random texture, or
//! Euler-disc-specific frequency is invented here.
//!
//! The footprint kernel is the exact one-dimensional marginal of the Hertz
//! elliptic pressure field along the declared travel direction:
//!
//! `w(x) = 3/(4a) * (1 - (x/a)^2), |x| <= a`,
//!
//! where `a` is the projected half-width of the ellipse. Uniform trace
//! samples are interpreted as a piecewise-linear measured profile, and both
//! the filtered height and its path derivative are integrated analytically on
//! every intersected segment. The only contact-law approximation in this leaf
//! is the explicit tangent linearization `delta_F = k_n * delta_h`.

use core::fmt;

use crate::{InputAuthority, InterfaceSystemRef};

const MAX_TRACE_SAMPLES: usize = 1_000_000;
const MAX_INTEGRATION_SEGMENTS: usize = 1_000_000;
const MAX_PERIODS_PER_FOOTPRINT: f64 = 8.0;

/// Boundary semantics of a uniformly sampled surface-height trace.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SurfaceTraceBoundary {
    /// Samples cover one finite track from `0` through `(n-1) * spacing`.
    Finite,
    /// Samples cover one period on `[0, n * spacing)`; the final segment wraps.
    Periodic,
}

/// Named, uniformly sampled, piecewise-linear surface height along one track.
///
/// Heights are positive out of the owning bulk surface. The texture-frame
/// identity must therefore also define the track direction and height sign.
#[derive(Clone, Debug, PartialEq)]
pub struct UniformSurfaceTrace {
    texture_frame_id: String,
    source_id: String,
    authority: InputAuthority,
    sample_spacing_m: f64,
    heights_m: Vec<f64>,
    boundary: SurfaceTraceBoundary,
}

impl UniformSurfaceTrace {
    /// Admit one measured or explicitly declared trace without upgrading its authority.
    pub fn new(
        texture_frame_id: impl Into<String>,
        source_id: impl Into<String>,
        authority: InputAuthority,
        sample_spacing_m: f64,
        heights_m: Vec<f64>,
        boundary: SurfaceTraceBoundary,
    ) -> Result<Self, SurfaceExcitationError> {
        let trace = Self {
            texture_frame_id: texture_frame_id.into(),
            source_id: source_id.into(),
            authority,
            sample_spacing_m,
            heights_m,
            boundary,
        };
        trace.validate()?;
        Ok(trace)
    }

    /// Opaque frame identity that defines the track and outward-height sign.
    #[must_use]
    pub fn texture_frame_id(&self) -> &str {
        &self.texture_frame_id
    }

    /// Caller source identity; this leaf does not validate the source.
    #[must_use]
    pub fn source_id(&self) -> &str {
        &self.source_id
    }

    /// Authority ceiling inherited from the caller.
    #[must_use]
    pub const fn authority(&self) -> InputAuthority {
        self.authority
    }

    /// Uniform spacing between adjacent retained samples [m].
    #[must_use]
    pub const fn sample_spacing_m(&self) -> f64 {
        self.sample_spacing_m
    }

    /// Retained outward-positive height samples [m].
    #[must_use]
    pub fn heights_m(&self) -> &[f64] {
        &self.heights_m
    }

    /// Declared finite or periodic boundary semantics.
    #[must_use]
    pub const fn boundary(&self) -> SurfaceTraceBoundary {
        self.boundary
    }

    /// Finite track length or periodic circumference [m].
    #[must_use]
    pub fn track_length_m(&self) -> f64 {
        let intervals = match self.boundary {
            SurfaceTraceBoundary::Finite => self.heights_m.len() - 1,
            SurfaceTraceBoundary::Periodic => self.heights_m.len(),
        };
        intervals as f64 * self.sample_spacing_m
    }

    fn validate(&self) -> Result<(), SurfaceExcitationError> {
        require_identity(&self.texture_frame_id, "texture_frame_id")?;
        require_identity(&self.source_id, "source_id")?;
        if !(self.sample_spacing_m.is_finite() && self.sample_spacing_m > 0.0) {
            return Err(SurfaceExcitationError::InvalidInput {
                field: "sample_spacing_m",
            });
        }
        let minimum = match self.boundary {
            SurfaceTraceBoundary::Finite => 2,
            SurfaceTraceBoundary::Periodic => 3,
        };
        if self.heights_m.len() < minimum || self.heights_m.len() > MAX_TRACE_SAMPLES {
            return Err(SurfaceExcitationError::TraceSampleCount {
                observed: self.heights_m.len(),
                minimum,
                maximum: MAX_TRACE_SAMPLES,
            });
        }
        if self.heights_m.iter().any(|height| !height.is_finite()) {
            return Err(SurfaceExcitationError::InvalidInput { field: "heights_m" });
        }
        let track_length_m = self.track_length_m();
        if !(track_length_m.is_finite() && track_length_m > 0.0) {
            return Err(SurfaceExcitationError::InvalidInput {
                field: "track_length_m",
            });
        }
        Ok(())
    }
}

/// Elliptic Hertz footprint and travel direction in its principal frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ProjectedHertzFootprint {
    /// Semimajor pressure-patch axis [m].
    pub semi_major_axis_m: f64,
    /// Semiminor pressure-patch axis [m].
    pub semi_minor_axis_m: f64,
    /// Travel direction measured from the major axis [rad].
    pub travel_angle_from_major_rad: f64,
}

impl ProjectedHertzFootprint {
    /// Exact support half-width of the elliptic footprint projected onto travel.
    pub fn projected_half_width_m(self) -> Result<f64, SurfaceExcitationError> {
        for (field, value) in [
            ("semi_major_axis_m", self.semi_major_axis_m),
            ("semi_minor_axis_m", self.semi_minor_axis_m),
        ] {
            if !(value.is_finite() && value > 0.0) {
                return Err(SurfaceExcitationError::InvalidInput { field });
            }
        }
        if !self.travel_angle_from_major_rad.is_finite() {
            return Err(SurfaceExcitationError::InvalidInput {
                field: "travel_angle_from_major_rad",
            });
        }
        let (sine, cosine) = self.travel_angle_from_major_rad.sin_cos();
        let half_width_m = (self.semi_major_axis_m * self.semi_major_axis_m * cosine * cosine
            + self.semi_minor_axis_m * self.semi_minor_axis_m * sine * sine)
            .sqrt();
        if !(half_width_m.is_finite() && half_width_m > 0.0) {
            return Err(SurfaceExcitationError::InvalidInput {
                field: "projected_half_width_m",
            });
        }
        Ok(half_width_m)
    }
}

/// One surface trace evaluated at its material-frame contact coordinate.
#[derive(Clone, Copy, Debug)]
pub struct SurfaceTraceMotion<'a> {
    /// Named height trace in the surface material frame.
    pub trace: &'a UniformSurfaceTrace,
    /// Current contact-path coordinate [m]. Periodic traces wrap exactly.
    pub path_coordinate_m: f64,
    /// Coordinate rate in that same texture frame [m/s].
    pub path_speed_m_per_s: f64,
}

/// Inputs for one tangent-linearized finite-patch roughness excitation sample.
#[derive(Clone, Copy, Debug)]
pub struct HertzRoughnessExcitationInput<'a> {
    /// Ordered, history-bearing dry interface identity.
    pub interface: &'a InterfaceSystemRef,
    /// First ordered surface.
    pub surface_a: SurfaceTraceMotion<'a>,
    /// Second ordered surface.
    pub surface_b: SurfaceTraceMotion<'a>,
    /// Resolved Hertz pressure footprint.
    pub footprint: ProjectedHertzFootprint,
    /// Nominal smooth-contact approach about which the law is linearized [m].
    pub nominal_approach_m: f64,
    /// Smooth-contact normal force at the linearization state [N].
    pub nominal_normal_force_n: f64,
    /// Consistent local derivative `dF_n / d approach` [N/m].
    pub normal_tangent_n_per_m: f64,
    /// Maximum admitted `abs(filtered_height) / nominal_approach`.
    pub maximum_linearized_height_fraction: f64,
}

/// Exact finite-patch filtering plus the explicitly linearized force result.
#[derive(Clone, Debug, PartialEq)]
pub struct HertzRoughnessExcitationReceipt {
    /// Ordered interface-system identity supplied to this leaf.
    pub ordered_interface_system_id: String,
    /// Texture-frame identities, in interface order.
    pub texture_frame_ids: [String; 2],
    /// Source identities, in interface order.
    pub source_ids: [String; 2],
    /// Caller authority ceilings, in interface order.
    pub authorities: [InputAuthority; 2],
    /// Projected Hertz pressure half-width along travel [m].
    pub projected_half_width_m: f64,
    /// Patch-filtered outward height for each surface [m].
    pub filtered_surface_heights_m: [f64; 2],
    /// Spatial derivative of each filtered height in its texture frame [m/m].
    pub filtered_surface_slopes: [f64; 2],
    /// Sum of both outward-positive filtered heights [m].
    pub combined_effective_height_m: f64,
    /// Material derivative of combined effective height [m/s].
    pub combined_effective_height_rate_m_per_s: f64,
    /// Tangent-linearized perturbation about the nominal normal load [N].
    pub normal_force_perturbation_n: f64,
    /// Time derivative of that perturbation [N/s].
    pub normal_force_perturbation_rate_n_per_s: f64,
    /// Actual admitted height/approach ratio.
    pub linearized_height_fraction: f64,
}

/// Evaluate deterministic surface excitation without synthesizing missing texture.
pub fn evaluate_hertz_roughness_excitation(
    input: HertzRoughnessExcitationInput<'_>,
) -> Result<HertzRoughnessExcitationReceipt, SurfaceExcitationError> {
    for (field, value) in [
        ("nominal_approach_m", input.nominal_approach_m),
        ("normal_tangent_n_per_m", input.normal_tangent_n_per_m),
    ] {
        if !(value.is_finite() && value > 0.0) {
            return Err(SurfaceExcitationError::InvalidInput { field });
        }
    }
    if !(input.nominal_normal_force_n.is_finite() && input.nominal_normal_force_n >= 0.0) {
        return Err(SurfaceExcitationError::InvalidInput {
            field: "nominal_normal_force_n",
        });
    }
    if !(input.maximum_linearized_height_fraction.is_finite()
        && input.maximum_linearized_height_fraction > 0.0
        && input.maximum_linearized_height_fraction <= 1.0)
    {
        return Err(SurfaceExcitationError::InvalidInput {
            field: "maximum_linearized_height_fraction",
        });
    }
    let half_width_m = input.footprint.projected_half_width_m()?;
    let a = filter_trace(input.surface_a, half_width_m)?;
    let b = filter_trace(input.surface_b, half_width_m)?;
    let combined_height_m = finite(a.height_m + b.height_m, "combined_effective_height_m")?;
    let combined_rate_m_per_s = finite(
        a.slope * input.surface_a.path_speed_m_per_s + b.slope * input.surface_b.path_speed_m_per_s,
        "combined_effective_height_rate_m_per_s",
    )?;
    let fraction = finite(
        combined_height_m.abs() / input.nominal_approach_m,
        "linearized_height_fraction",
    )?;
    if fraction > input.maximum_linearized_height_fraction {
        return Err(SurfaceExcitationError::OutsideLinearizedContact {
            observed_fraction: fraction,
            maximum_fraction: input.maximum_linearized_height_fraction,
        });
    }
    let force_perturbation_n = finite(
        input.normal_tangent_n_per_m * combined_height_m,
        "normal_force_perturbation_n",
    )?;
    let perturbed_force_n = finite(
        input.nominal_normal_force_n + force_perturbation_n,
        "perturbed_normal_force_n",
    )?;
    if perturbed_force_n < 0.0 {
        return Err(SurfaceExcitationError::WouldOpenContact { perturbed_force_n });
    }
    let force_rate_n_per_s = finite(
        input.normal_tangent_n_per_m * combined_rate_m_per_s,
        "normal_force_perturbation_rate_n_per_s",
    )?;
    Ok(HertzRoughnessExcitationReceipt {
        ordered_interface_system_id: input.interface.ordered_system_id().to_owned(),
        texture_frame_ids: [
            input.surface_a.trace.texture_frame_id.clone(),
            input.surface_b.trace.texture_frame_id.clone(),
        ],
        source_ids: [
            input.surface_a.trace.source_id.clone(),
            input.surface_b.trace.source_id.clone(),
        ],
        authorities: [
            input.surface_a.trace.authority,
            input.surface_b.trace.authority,
        ],
        projected_half_width_m: half_width_m,
        filtered_surface_heights_m: [a.height_m, b.height_m],
        filtered_surface_slopes: [a.slope, b.slope],
        combined_effective_height_m: combined_height_m,
        combined_effective_height_rate_m_per_s: combined_rate_m_per_s,
        normal_force_perturbation_n: force_perturbation_n,
        normal_force_perturbation_rate_n_per_s: force_rate_n_per_s,
        linearized_height_fraction: fraction,
    })
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct FilteredTrace {
    height_m: f64,
    slope: f64,
}

fn filter_trace(
    motion: SurfaceTraceMotion<'_>,
    half_width_m: f64,
) -> Result<FilteredTrace, SurfaceExcitationError> {
    if !(motion.path_coordinate_m.is_finite() && motion.path_speed_m_per_s.is_finite()) {
        return Err(SurfaceExcitationError::InvalidInput {
            field: "surface path coordinate or speed",
        });
    }
    motion.trace.validate()?;
    let track_length_m = motion.trace.track_length_m();
    let center_m = match motion.trace.boundary {
        SurfaceTraceBoundary::Finite => {
            let tolerance = 64.0 * f64::EPSILON * track_length_m.max(1.0);
            if motion.path_coordinate_m - half_width_m < -tolerance
                || motion.path_coordinate_m + half_width_m > track_length_m + tolerance
            {
                return Err(SurfaceExcitationError::FootprintOutsideFiniteTrace {
                    center_m: motion.path_coordinate_m,
                    half_width_m,
                    track_length_m,
                });
            }
            motion.path_coordinate_m
        }
        SurfaceTraceBoundary::Periodic => {
            if 2.0 * half_width_m > MAX_PERIODS_PER_FOOTPRINT * track_length_m {
                return Err(SurfaceExcitationError::FootprintTooWideForPeriodicTrace {
                    footprint_width_m: 2.0 * half_width_m,
                    track_length_m,
                });
            }
            motion.path_coordinate_m.rem_euclid(track_length_m)
        }
    };
    let spacing_m = motion.trace.sample_spacing_m;
    let (start_m, end_m) = match motion.trace.boundary {
        SurfaceTraceBoundary::Finite => (
            (center_m - half_width_m).max(0.0),
            (center_m + half_width_m).min(track_length_m),
        ),
        SurfaceTraceBoundary::Periodic => (center_m - half_width_m, center_m + half_width_m),
    };
    let first_segment = (start_m / spacing_m).floor() as i64;
    let last_segment = (end_m / spacing_m).ceil() as i64 - 1;
    let segment_count_i64 = last_segment.saturating_sub(first_segment).saturating_add(1);
    let segment_count = usize::try_from(segment_count_i64).map_err(|_| {
        SurfaceExcitationError::IntegrationBudgetExceeded {
            requested_segments: usize::MAX,
            maximum_segments: MAX_INTEGRATION_SEGMENTS,
        }
    })?;
    if segment_count > MAX_INTEGRATION_SEGMENTS {
        return Err(SurfaceExcitationError::IntegrationBudgetExceeded {
            requested_segments: segment_count,
            maximum_segments: MAX_INTEGRATION_SEGMENTS,
        });
    }

    let mut height_m = 0.0;
    let mut slope = 0.0;
    for segment in first_segment..=last_segment {
        let segment_start_m = segment as f64 * spacing_m;
        let lo_m = start_m.max(segment_start_m);
        let hi_m = end_m.min(segment_start_m + spacing_m);
        if hi_m <= lo_m {
            continue;
        }
        let (height_start_m, height_end_m) = segment_heights(motion.trace, segment)?;
        let profile_slope = (height_end_m - height_start_m) / spacing_m;
        let height_at_center_m = height_start_m + profile_slope * (center_m - segment_start_m);
        let u0 = lo_m - center_m;
        let u1 = hi_m - center_m;
        let h2 = half_width_m * half_width_m;
        let primitive_height = |u: f64| {
            height_at_center_m * (u - u.powi(3) / (3.0 * h2))
                + profile_slope * (0.5 * u * u - u.powi(4) / (4.0 * h2))
        };
        let primitive_slope =
            |u: f64| height_at_center_m * 0.5 * u * u + profile_slope * u.powi(3) / 3.0;
        height_m += 3.0 / (4.0 * half_width_m) * (primitive_height(u1) - primitive_height(u0));
        slope += 3.0 / (2.0 * half_width_m.powi(3)) * (primitive_slope(u1) - primitive_slope(u0));
    }
    Ok(FilteredTrace {
        height_m: finite(height_m, "filtered_surface_height_m")?,
        slope: finite(slope, "filtered_surface_slope")?,
    })
}

fn segment_heights(
    trace: &UniformSurfaceTrace,
    segment: i64,
) -> Result<(f64, f64), SurfaceExcitationError> {
    match trace.boundary {
        SurfaceTraceBoundary::Finite => {
            let start = usize::try_from(segment).map_err(|_| {
                SurfaceExcitationError::FootprintOutsideFiniteTrace {
                    center_m: 0.0,
                    half_width_m: 0.0,
                    track_length_m: trace.track_length_m(),
                }
            })?;
            let end = start
                .checked_add(1)
                .ok_or(SurfaceExcitationError::InvalidInput {
                    field: "finite trace segment",
                })?;
            let Some((&a, &b)) = trace.heights_m.get(start).zip(trace.heights_m.get(end)) else {
                return Err(SurfaceExcitationError::FootprintOutsideFiniteTrace {
                    center_m: 0.0,
                    half_width_m: 0.0,
                    track_length_m: trace.track_length_m(),
                });
            };
            Ok((a, b))
        }
        SurfaceTraceBoundary::Periodic => {
            let count = i64::try_from(trace.heights_m.len()).map_err(|_| {
                SurfaceExcitationError::InvalidInput {
                    field: "periodic trace sample count",
                }
            })?;
            let start = usize::try_from(segment.rem_euclid(count)).map_err(|_| {
                SurfaceExcitationError::InvalidInput {
                    field: "periodic trace segment",
                }
            })?;
            let end = (start + 1) % trace.heights_m.len();
            Ok((trace.heights_m[start], trace.heights_m[end]))
        }
    }
}

fn require_identity(value: &str, field: &'static str) -> Result<(), SurfaceExcitationError> {
    if value.trim().is_empty() || value.len() > 1024 {
        return Err(SurfaceExcitationError::MissingIdentity { field });
    }
    Ok(())
}

fn finite(value: f64, field: &'static str) -> Result<f64, SurfaceExcitationError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SurfaceExcitationError::InvalidInput { field })
    }
}

/// Total refusal surface for contact-patch-filtered surface excitation.
#[derive(Clone, Debug, PartialEq)]
pub enum SurfaceExcitationError {
    /// A texture/source identity was blank or unreasonably large.
    MissingIdentity {
        /// Stable offending field name.
        field: &'static str,
    },
    /// A physical scalar, coordinate, or derived value was invalid.
    InvalidInput {
        /// Stable offending field name.
        field: &'static str,
    },
    /// The trace was too short for its boundary semantics or exceeded its cap.
    TraceSampleCount {
        /// Supplied sample count.
        observed: usize,
        /// Minimum admitted by the selected boundary semantics.
        minimum: usize,
        /// Resource ceiling.
        maximum: usize,
    },
    /// A finite trace does not cover the complete Hertz footprint.
    FootprintOutsideFiniteTrace {
        /// Offered footprint center [m].
        center_m: f64,
        /// Offered footprint half-width [m].
        half_width_m: f64,
        /// Available finite track length [m].
        track_length_m: f64,
    },
    /// A periodic footprint spans too many repeated copies for this scalar leaf.
    FootprintTooWideForPeriodicTrace {
        /// Complete projected footprint width [m].
        footprint_width_m: f64,
        /// Declared spatial period [m].
        track_length_m: f64,
    },
    /// The exact piecewise integration would exceed its bounded work cap.
    IntegrationBudgetExceeded {
        /// Number of piecewise-linear segments required.
        requested_segments: usize,
        /// Fixed scalar-leaf work ceiling.
        maximum_segments: usize,
    },
    /// The declared contact-tangent linearization is being used too far from its state.
    OutsideLinearizedContact {
        /// Actual absolute filtered-height/approach ratio.
        observed_fraction: f64,
        /// Caller-declared admissible ceiling.
        maximum_fraction: f64,
    },
    /// The linearized perturbation would make the nominal unilateral load negative.
    WouldOpenContact {
        /// Nominal plus perturbation normal force [N].
        perturbed_force_n: f64,
    },
}

impl fmt::Display for SurfaceExcitationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for SurfaceExcitationError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::InterfaceMedium;

    fn interface() -> InterfaceSystemRef {
        InterfaceSystemRef::new(
            "disc-edge->glass",
            "measured-run-in-state",
            "interface-card-hash",
            InputAuthority::CallerDeclared,
            InterfaceMedium::Dry,
        )
        .expect("valid dry interface")
    }

    fn trace(
        id: &str,
        spacing_m: f64,
        heights_m: Vec<f64>,
        boundary: SurfaceTraceBoundary,
    ) -> UniformSurfaceTrace {
        UniformSurfaceTrace::new(
            id,
            format!("source/{id}"),
            InputAuthority::SyntheticFixture,
            spacing_m,
            heights_m,
            boundary,
        )
        .expect("valid trace")
    }

    #[test]
    fn g1_affine_trace_is_filtered_and_differentiated_exactly() {
        let spacing_m = 0.001;
        let slope = 2.5e-4;
        let intercept_m = 3.0e-6;
        let samples = (0..101)
            .map(|index| intercept_m + slope * index as f64 * spacing_m)
            .collect();
        let affine = trace("affine", spacing_m, samples, SurfaceTraceBoundary::Finite);
        let flat = trace(
            "flat",
            spacing_m,
            vec![0.0; 101],
            SurfaceTraceBoundary::Finite,
        );
        let coordinate_m = 0.0473;
        let interface = interface();
        let receipt = evaluate_hertz_roughness_excitation(HertzRoughnessExcitationInput {
            interface: &interface,
            surface_a: SurfaceTraceMotion {
                trace: &affine,
                path_coordinate_m: coordinate_m,
                path_speed_m_per_s: 0.7,
            },
            surface_b: SurfaceTraceMotion {
                trace: &flat,
                path_coordinate_m: coordinate_m,
                path_speed_m_per_s: 0.0,
            },
            footprint: ProjectedHertzFootprint {
                semi_major_axis_m: 0.003,
                semi_minor_axis_m: 0.001,
                travel_angle_from_major_rad: 0.0,
            },
            nominal_approach_m: 1.0e-3,
            nominal_normal_force_n: 100.0,
            normal_tangent_n_per_m: 2.0e6,
            maximum_linearized_height_fraction: 0.1,
        })
        .expect("affine trace lies in the linearized domain");
        let expected_height_m = intercept_m + slope * coordinate_m;
        assert!((receipt.filtered_surface_heights_m[0] - expected_height_m).abs() < 2.0e-18);
        assert!((receipt.filtered_surface_slopes[0] - slope).abs() < 2.0e-15);
        assert!((receipt.combined_effective_height_rate_m_per_s - slope * 0.7).abs() < 2.0e-15);
        assert!((receipt.normal_force_perturbation_n - 2.0e6 * expected_height_m).abs() < 2.0e-11);
    }

    #[test]
    fn g0_constant_periodic_trace_wraps_without_a_seam_or_force_rate() {
        let constant = trace(
            "constant-periodic",
            0.001,
            vec![2.0e-6; 32],
            SurfaceTraceBoundary::Periodic,
        );
        let zero = trace(
            "zero-periodic",
            0.001,
            vec![0.0; 32],
            SurfaceTraceBoundary::Periodic,
        );
        let interface = interface();
        let make = |coordinate_m| HertzRoughnessExcitationInput {
            interface: &interface,
            surface_a: SurfaceTraceMotion {
                trace: &constant,
                path_coordinate_m: coordinate_m,
                path_speed_m_per_s: 1.0,
            },
            surface_b: SurfaceTraceMotion {
                trace: &zero,
                path_coordinate_m: 0.0,
                path_speed_m_per_s: 0.0,
            },
            footprint: ProjectedHertzFootprint {
                semi_major_axis_m: 0.002,
                semi_minor_axis_m: 0.001,
                travel_angle_from_major_rad: core::f64::consts::FRAC_PI_2,
            },
            nominal_approach_m: 1.0e-3,
            nominal_normal_force_n: 100.0,
            normal_tangent_n_per_m: 1.0e6,
            maximum_linearized_height_fraction: 0.1,
        };
        let before = evaluate_hertz_roughness_excitation(make(0.0319)).expect("before wrap");
        let after = evaluate_hertz_roughness_excitation(make(0.0639)).expect("period translated");
        assert!((before.combined_effective_height_m - 2.0e-6).abs() < 2.0e-19);
        assert!(before.combined_effective_height_rate_m_per_s.abs() < 2.0e-15);
        assert_eq!(
            before.filtered_surface_heights_m,
            after.filtered_surface_heights_m
        );
        assert_eq!(
            before.filtered_surface_slopes,
            after.filtered_surface_slopes
        );
    }

    #[test]
    fn g3_projected_ellipse_axes_swap_under_quarter_turn() {
        let major = ProjectedHertzFootprint {
            semi_major_axis_m: 0.003,
            semi_minor_axis_m: 0.001,
            travel_angle_from_major_rad: 0.0,
        }
        .projected_half_width_m()
        .expect("major projection");
        let minor = ProjectedHertzFootprint {
            semi_major_axis_m: 0.003,
            semi_minor_axis_m: 0.001,
            travel_angle_from_major_rad: core::f64::consts::FRAC_PI_2,
        }
        .projected_half_width_m()
        .expect("minor projection");
        assert!((major - 0.003).abs() < 1.0e-18);
        assert!((minor - 0.001).abs() < 1.0e-18);
    }

    #[test]
    fn g0_missing_texture_and_out_of_domain_linearization_refuse() {
        assert!(matches!(
            UniformSurfaceTrace::new(
                " ",
                "source",
                InputAuthority::CallerDeclared,
                1.0,
                vec![0.0, 0.0],
                SurfaceTraceBoundary::Finite,
            ),
            Err(SurfaceExcitationError::MissingIdentity {
                field: "texture_frame_id"
            })
        ));

        let high = trace(
            "high",
            0.001,
            vec![2.0e-4; 32],
            SurfaceTraceBoundary::Periodic,
        );
        let zero = trace("zero", 0.001, vec![0.0; 32], SurfaceTraceBoundary::Periodic);
        let interface = interface();
        let result = evaluate_hertz_roughness_excitation(HertzRoughnessExcitationInput {
            interface: &interface,
            surface_a: SurfaceTraceMotion {
                trace: &high,
                path_coordinate_m: 0.0,
                path_speed_m_per_s: 0.0,
            },
            surface_b: SurfaceTraceMotion {
                trace: &zero,
                path_coordinate_m: 0.0,
                path_speed_m_per_s: 0.0,
            },
            footprint: ProjectedHertzFootprint {
                semi_major_axis_m: 0.001,
                semi_minor_axis_m: 0.001,
                travel_angle_from_major_rad: 0.0,
            },
            nominal_approach_m: 1.0e-3,
            nominal_normal_force_n: 100.0,
            normal_tangent_n_per_m: 1.0e6,
            maximum_linearized_height_fraction: 0.1,
        });
        assert!(matches!(
            result,
            Err(SurfaceExcitationError::OutsideLinearizedContact { .. })
        ));
    }
}

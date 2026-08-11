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
use std::sync::Arc;

use fs_math::det;
use fs_rand::StreamKey;

use crate::{InputAuthority, InterfaceSystemRef};

const MAX_TRACE_SAMPLES: usize = 1_000_000;
const MAX_INTEGRATION_SEGMENTS: usize = 1_000_000;
const MAX_PERIODS_PER_FOOTPRINT: f64 = 8.0;
const MIN_SAMPLES_PER_SHORTEST_HARMONIC: usize = 8;
const SELF_AFFINE_PHASE_KERNEL: u32 = 0x5341_4646;

/// One real Fourier component of a periodic measured or declared surface track.
///
/// `cosine_amplitude_m` and `sine_amplitude_m` are signed outward-height
/// coefficients. The integer cycle count makes the generated trace exactly
/// periodic on its declared material-frame track; no windowing or hidden seam
/// correction is applied.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PeriodicSurfaceHarmonic {
    /// Positive integer cycles around the complete material-frame track.
    pub cycles_per_track: u32,
    /// Cosine coefficient [m].
    pub cosine_amplitude_m: f64,
    /// Sine coefficient [m].
    pub sine_amplitude_m: f64,
}

/// Bounded random-phase realization of a one-dimensional self-affine profile.
///
/// The retained spatial power follows `Phi(q) proportional to q^-(2 H + 1)`
/// between the declared integer cycle cutoffs. Fourier magnitudes are
/// normalized so the continuous zero-mean periodic profile has exactly the
/// requested RMS height (up to floating-point roundoff). The seed chooses only
/// spatial phases; it never selects a temporal or acoustic frequency.
///
/// This is a caller-parameterized geometry model. It is not a measurement, a
/// material-name preset, a two-dimensional areal texture, or a statistical
/// validation of any specimen.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SelfAffinePeriodicProfileSpectrum {
    rms_height_m: f64,
    hurst_exponent: f64,
    minimum_cycles_per_track: u32,
    maximum_cycles_per_track: u32,
    phase_seed: u64,
}

impl SelfAffinePeriodicProfileSpectrum {
    /// Admit one band-limited self-affine profile spectrum.
    ///
    /// `0 < hurst_exponent < 1` is the ordinary self-affine profile range.
    /// Both cycle cutoffs are positive and inclusive. The returned spectrum is
    /// still only caller input; construction does not upgrade its authority.
    pub fn new(
        rms_height_m: f64,
        hurst_exponent: f64,
        minimum_cycles_per_track: u32,
        maximum_cycles_per_track: u32,
        phase_seed: u64,
    ) -> Result<Self, SurfaceExcitationError> {
        if !(rms_height_m.is_finite() && rms_height_m > 0.0) {
            return Err(SurfaceExcitationError::InvalidInput {
                field: "self_affine_rms_height_m",
            });
        }
        if !(hurst_exponent.is_finite() && hurst_exponent > 0.0 && hurst_exponent < 1.0) {
            return Err(SurfaceExcitationError::InvalidInput {
                field: "self_affine_hurst_exponent",
            });
        }
        if minimum_cycles_per_track == 0 || maximum_cycles_per_track < minimum_cycles_per_track {
            return Err(SurfaceExcitationError::InvalidInput {
                field: "self_affine_cycle_band",
            });
        }
        let maximum_resolvable_cycle = MAX_TRACE_SAMPLES / MIN_SAMPLES_PER_SHORTEST_HARMONIC;
        let maximum_cycle = usize::try_from(maximum_cycles_per_track).map_err(|_| {
            SurfaceExcitationError::SurfaceSpectrumUnderresolved {
                sample_count: MAX_TRACE_SAMPLES,
                required_samples: usize::MAX,
            }
        })?;
        if maximum_cycle > maximum_resolvable_cycle {
            return Err(SurfaceExcitationError::SurfaceSpectrumUnderresolved {
                sample_count: MAX_TRACE_SAMPLES,
                required_samples: maximum_cycle
                    .checked_mul(MIN_SAMPLES_PER_SHORTEST_HARMONIC)
                    .unwrap_or(usize::MAX),
            });
        }
        let harmonic_count = u64::from(maximum_cycles_per_track)
            .checked_sub(u64::from(minimum_cycles_per_track))
            .and_then(|span| span.checked_add(1))
            .and_then(|count| usize::try_from(count).ok())
            .ok_or(SurfaceExcitationError::TraceCapacity {
                requested: usize::MAX,
            })?;
        if harmonic_count > MAX_TRACE_SAMPLES {
            return Err(SurfaceExcitationError::TraceCapacity {
                requested: harmonic_count,
            });
        }
        Ok(Self {
            rms_height_m,
            hurst_exponent,
            minimum_cycles_per_track,
            maximum_cycles_per_track,
            phase_seed,
        })
    }

    /// Requested continuous-profile RMS height [m].
    #[must_use]
    pub const fn rms_height_m(&self) -> f64 {
        self.rms_height_m
    }

    /// One-dimensional self-affine Hurst exponent.
    #[must_use]
    pub const fn hurst_exponent(&self) -> f64 {
        self.hurst_exponent
    }

    /// Inclusive lowest retained spatial cycle around the track.
    #[must_use]
    pub const fn minimum_cycles_per_track(&self) -> u32 {
        self.minimum_cycles_per_track
    }

    /// Inclusive highest retained spatial cycle around the track.
    #[must_use]
    pub const fn maximum_cycles_per_track(&self) -> u32 {
        self.maximum_cycles_per_track
    }

    /// Explicit deterministic phase-realization seed.
    #[must_use]
    pub const fn phase_seed(&self) -> u64 {
        self.phase_seed
    }

    /// Realize the declared PSD into explicit periodic Fourier coefficients.
    ///
    /// Each cycle owns a separate counter-based random stream, so extending or
    /// truncating the admitted band does not change the phases of shared
    /// cycles. Coefficient order is canonical ascending cycle order.
    pub fn realize_harmonics(
        &self,
    ) -> Result<Vec<PeriodicSurfaceHarmonic>, SurfaceExcitationError> {
        let count = usize::try_from(
            u64::from(self.maximum_cycles_per_track) - u64::from(self.minimum_cycles_per_track) + 1,
        )
        .map_err(|_| SurfaceExcitationError::TraceCapacity {
            requested: usize::MAX,
        })?;
        let spectral_exponent = self.hurst_exponent.mul_add(2.0, 1.0);
        let mut weights = Vec::new();
        weights
            .try_reserve_exact(count)
            .map_err(|_| SurfaceExcitationError::TraceCapacity { requested: count })?;
        let mut total_weight = 0.0_f64;
        for cycles in self.minimum_cycles_per_track..=self.maximum_cycles_per_track {
            let cycle = f64::from(cycles);
            let weight = det::exp(-spectral_exponent * det::ln(cycle));
            if !(weight.is_finite() && weight > 0.0) {
                return Err(SurfaceExcitationError::InvalidInput {
                    field: "self_affine_spectral_weight",
                });
            }
            total_weight += weight;
            weights.push((cycles, weight));
        }
        if !(total_weight.is_finite() && total_weight > 0.0) {
            return Err(SurfaceExcitationError::InvalidInput {
                field: "self_affine_total_spectral_weight",
            });
        }

        let mut harmonics = Vec::new();
        harmonics
            .try_reserve_exact(count)
            .map_err(|_| SurfaceExcitationError::TraceCapacity { requested: count })?;
        for (cycles, weight) in weights {
            let amplitude_m = self.rms_height_m * det::sqrt(2.0 * weight / total_weight);
            let mut phase_stream = StreamKey {
                seed: self.phase_seed,
                kernel: SELF_AFFINE_PHASE_KERNEL,
                tile: cycles,
            }
            .stream();
            let phase_rad = core::f64::consts::TAU * phase_stream.next_f64();
            harmonics.push(PeriodicSurfaceHarmonic {
                cycles_per_track: cycles,
                cosine_amplitude_m: amplitude_m * det::cos(phase_rad),
                sine_amplitude_m: amplitude_m * det::sin(phase_rad),
            });
        }
        Ok(harmonics)
    }
}

/// Explicit spatial spectrum from which a periodic surface-height trace is sampled.
///
/// This is a geometry input, not a sound synthesizer or a material-name preset.
/// A measured Fourier decomposition can be supplied directly; a hypothetical
/// surface may use caller-declared coefficients under an appropriately weak
/// authority. Components are canonicalized by cycle count, so their input
/// order cannot change the sampled trace.
#[derive(Clone, Debug, PartialEq)]
pub struct PeriodicHarmonicSurface {
    texture_frame_id: String,
    source_id: String,
    authority: InputAuthority,
    track_length_m: f64,
    sample_count: usize,
    harmonics: Vec<PeriodicSurfaceHarmonic>,
}

impl PeriodicHarmonicSurface {
    /// Admit an explicit periodic surface spectrum.
    ///
    /// The sample grid must retain at least eight samples per period of the
    /// shortest harmonic. Duplicate cycle counts refuse rather than depending
    /// on caller summation order. An empty harmonic set is a valid declared
    /// perfectly smooth track and realizes to exact zero heights.
    pub fn new(
        texture_frame_id: impl Into<String>,
        source_id: impl Into<String>,
        authority: InputAuthority,
        track_length_m: f64,
        sample_count: usize,
        mut harmonics: Vec<PeriodicSurfaceHarmonic>,
    ) -> Result<Self, SurfaceExcitationError> {
        let texture_frame_id = texture_frame_id.into();
        let source_id = source_id.into();
        require_identity(&texture_frame_id, "texture_frame_id")?;
        require_identity(&source_id, "source_id")?;
        if !(track_length_m.is_finite() && track_length_m > 0.0) {
            return Err(SurfaceExcitationError::InvalidInput {
                field: "track_length_m",
            });
        }
        if !(3..=MAX_TRACE_SAMPLES).contains(&sample_count) {
            return Err(SurfaceExcitationError::TraceSampleCount {
                observed: sample_count,
                minimum: 3,
                maximum: MAX_TRACE_SAMPLES,
            });
        }
        for harmonic in &harmonics {
            if harmonic.cycles_per_track == 0
                || !harmonic.cosine_amplitude_m.is_finite()
                || !harmonic.sine_amplitude_m.is_finite()
            {
                return Err(SurfaceExcitationError::InvalidInput {
                    field: "periodic_surface_harmonic",
                });
            }
        }
        harmonics.sort_by_key(|harmonic| harmonic.cycles_per_track);
        if harmonics
            .windows(2)
            .any(|pair| pair[0].cycles_per_track == pair[1].cycles_per_track)
        {
            return Err(SurfaceExcitationError::DuplicateHarmonicCycle);
        }
        if let Some(maximum_cycle) = harmonics.last().map(|harmonic| harmonic.cycles_per_track) {
            let required_samples = usize::try_from(maximum_cycle)
                .ok()
                .and_then(|cycles| cycles.checked_mul(MIN_SAMPLES_PER_SHORTEST_HARMONIC))
                .ok_or(SurfaceExcitationError::SurfaceSpectrumUnderresolved {
                    sample_count,
                    required_samples: usize::MAX,
                })?;
            if sample_count < required_samples {
                return Err(SurfaceExcitationError::SurfaceSpectrumUnderresolved {
                    sample_count,
                    required_samples,
                });
            }
        }
        Ok(Self {
            texture_frame_id,
            source_id,
            authority,
            track_length_m,
            sample_count,
            harmonics,
        })
    }

    /// Exact material-frame period length [m].
    #[must_use]
    pub const fn track_length_m(&self) -> f64 {
        self.track_length_m
    }

    /// Canonically ordered spatial Fourier components.
    #[must_use]
    pub fn harmonics(&self) -> &[PeriodicSurfaceHarmonic] {
        &self.harmonics
    }

    /// Sample the explicit spectrum into the trace consumed by the finite-patch filter.
    pub fn realize(&self) -> Result<UniformSurfaceTrace, SurfaceExcitationError> {
        let mut heights_m = Vec::new();
        heights_m
            .try_reserve_exact(self.sample_count)
            .map_err(|_| SurfaceExcitationError::TraceCapacity {
                requested: self.sample_count,
            })?;
        let sample_count = self.sample_count as f64;
        for sample in 0..self.sample_count {
            let phase_fraction = sample as f64 / sample_count;
            let mut height_m = 0.0_f64;
            for harmonic in &self.harmonics {
                let phase =
                    core::f64::consts::TAU * f64::from(harmonic.cycles_per_track) * phase_fraction;
                let (sine, cosine) = phase.sin_cos();
                height_m = harmonic.cosine_amplitude_m.mul_add(cosine, height_m);
                height_m = harmonic.sine_amplitude_m.mul_add(sine, height_m);
            }
            if !height_m.is_finite() {
                return Err(SurfaceExcitationError::InvalidInput {
                    field: "realized_surface_height_m",
                });
            }
            heights_m.push(height_m);
        }
        UniformSurfaceTrace::new(
            self.texture_frame_id.clone(),
            self.source_id.clone(),
            self.authority,
            self.track_length_m / sample_count,
            heights_m,
            SurfaceTraceBoundary::Periodic,
        )
    }
}

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
    heights_m: Arc<[f64]>,
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
            heights_m: heights_m.into(),
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

/// Contact-path height and material derivative for an ordered surface pair.
///
/// This is a geometry receipt, not a force or sound model. A zero projected
/// half-width is the exact point-contact limit used to decide first touch;
/// positive half-widths use the Hertz pressure marginal documented above.
#[derive(Clone, Debug, PartialEq)]
pub struct FilteredSurfacePairReceipt {
    /// Ordered interface-system identity supplied to this leaf.
    pub ordered_interface_system_id: String,
    /// Texture-frame identities, in interface order.
    pub texture_frame_ids: [String; 2],
    /// Source identities, in interface order.
    pub source_ids: [String; 2],
    /// Caller authority ceilings, in interface order.
    pub authorities: [InputAuthority; 2],
    /// Projected pressure-support half-width [m], or zero for the point limit.
    pub projected_half_width_m: f64,
    /// Outward-positive height of each surface after the declared filtering [m].
    pub filtered_surface_heights_m: [f64; 2],
    /// Material-frame spatial derivative of each filtered height [m/m].
    pub filtered_surface_slopes: [f64; 2],
    /// Sum of both outward-positive heights [m].
    pub combined_effective_height_m: f64,
    /// Material derivative of the combined height along both paths [m/s].
    pub combined_effective_height_rate_m_per_s: f64,
}

/// Evaluate the exact point-support limit of two declared surface traces.
///
/// Piecewise-linear trace interpolation supplies both height and slope. This
/// limit is needed by unilateral contact event selection before a finite Hertz
/// footprint exists; it does not invent an asperity radius or contact force.
pub fn evaluate_point_surface_pair(
    interface: &InterfaceSystemRef,
    surface_a: SurfaceTraceMotion<'_>,
    surface_b: SurfaceTraceMotion<'_>,
) -> Result<FilteredSurfacePairReceipt, SurfaceExcitationError> {
    let a = sample_trace(surface_a)?;
    let b = sample_trace(surface_b)?;
    surface_pair_receipt(interface, surface_a, surface_b, 0.0, a, b)
}

/// Filter two declared surface traces through a resolved Hertz footprint.
///
/// The result remains purely geometric. A caller may use it to re-resolve a
/// nonlinear unilateral contact law, or may separately request the explicitly
/// tangent-linear force approximation below.
pub fn evaluate_hertz_filtered_surface_pair(
    interface: &InterfaceSystemRef,
    surface_a: SurfaceTraceMotion<'_>,
    surface_b: SurfaceTraceMotion<'_>,
    footprint: ProjectedHertzFootprint,
) -> Result<FilteredSurfacePairReceipt, SurfaceExcitationError> {
    let half_width_m = footprint.projected_half_width_m()?;
    let a = filter_trace(surface_a, half_width_m)?;
    let b = filter_trace(surface_b, half_width_m)?;
    surface_pair_receipt(interface, surface_a, surface_b, half_width_m, a, b)
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
    let filtered = evaluate_hertz_filtered_surface_pair(
        input.interface,
        input.surface_a,
        input.surface_b,
        input.footprint,
    )?;
    let combined_height_m = filtered.combined_effective_height_m;
    let combined_rate_m_per_s = filtered.combined_effective_height_rate_m_per_s;
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
        ordered_interface_system_id: filtered.ordered_interface_system_id,
        texture_frame_ids: filtered.texture_frame_ids,
        source_ids: filtered.source_ids,
        authorities: filtered.authorities,
        projected_half_width_m: filtered.projected_half_width_m,
        filtered_surface_heights_m: filtered.filtered_surface_heights_m,
        filtered_surface_slopes: filtered.filtered_surface_slopes,
        combined_effective_height_m: combined_height_m,
        combined_effective_height_rate_m_per_s: combined_rate_m_per_s,
        normal_force_perturbation_n: force_perturbation_n,
        normal_force_perturbation_rate_n_per_s: force_rate_n_per_s,
        linearized_height_fraction: fraction,
    })
}

fn surface_pair_receipt(
    interface: &InterfaceSystemRef,
    surface_a: SurfaceTraceMotion<'_>,
    surface_b: SurfaceTraceMotion<'_>,
    projected_half_width_m: f64,
    a: FilteredTrace,
    b: FilteredTrace,
) -> Result<FilteredSurfacePairReceipt, SurfaceExcitationError> {
    let combined_effective_height_m =
        finite(a.height_m + b.height_m, "combined_effective_height_m")?;
    let combined_effective_height_rate_m_per_s = finite(
        a.slope * surface_a.path_speed_m_per_s + b.slope * surface_b.path_speed_m_per_s,
        "combined_effective_height_rate_m_per_s",
    )?;
    Ok(FilteredSurfacePairReceipt {
        ordered_interface_system_id: interface.ordered_system_id().to_owned(),
        texture_frame_ids: [
            surface_a.trace.texture_frame_id.clone(),
            surface_b.trace.texture_frame_id.clone(),
        ],
        source_ids: [
            surface_a.trace.source_id.clone(),
            surface_b.trace.source_id.clone(),
        ],
        authorities: [surface_a.trace.authority, surface_b.trace.authority],
        projected_half_width_m,
        filtered_surface_heights_m: [a.height_m, b.height_m],
        filtered_surface_slopes: [a.slope, b.slope],
        combined_effective_height_m,
        combined_effective_height_rate_m_per_s,
    })
}

fn sample_trace(motion: SurfaceTraceMotion<'_>) -> Result<FilteredTrace, SurfaceExcitationError> {
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
            if motion.path_coordinate_m < -tolerance
                || motion.path_coordinate_m > track_length_m + tolerance
            {
                return Err(SurfaceExcitationError::FootprintOutsideFiniteTrace {
                    center_m: motion.path_coordinate_m,
                    half_width_m: 0.0,
                    track_length_m,
                });
            }
            motion.path_coordinate_m.clamp(0.0, track_length_m)
        }
        SurfaceTraceBoundary::Periodic => motion.path_coordinate_m.rem_euclid(track_length_m),
    };
    let spacing_m = motion.trace.sample_spacing_m;
    let segment = match motion.trace.boundary {
        SurfaceTraceBoundary::Finite if center_m == track_length_m => {
            i64::try_from(motion.trace.heights_m.len() - 2).map_err(|_| {
                SurfaceExcitationError::InvalidInput {
                    field: "finite trace point segment",
                }
            })?
        }
        _ => (center_m / spacing_m).floor() as i64,
    };
    let segment_start_m = segment as f64 * spacing_m;
    let (height_start_m, height_end_m) = segment_heights(motion.trace, segment)?;
    let slope = finite(
        (height_end_m - height_start_m) / spacing_m,
        "point_surface_slope",
    )?;
    let height_m = finite(
        height_start_m + slope * (center_m - segment_start_m),
        "point_surface_height_m",
    )?;
    Ok(FilteredTrace { height_m, slope })
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
                + profile_slope * (0.5 * u * u - det::powi(u, 4) / (4.0 * h2))
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
    /// Realizing a declared spectrum could not reserve its bounded trace.
    TraceCapacity {
        /// Requested number of retained samples.
        requested: usize,
    },
    /// Two explicit components named the same spatial cycle.
    DuplicateHarmonicCycle,
    /// The requested trace grid would alias or poorly resolve its shortest period.
    SurfaceSpectrumUnderresolved {
        /// Caller-selected periodic sample count.
        sample_count: usize,
        /// Minimum eight-samples-per-period count for the highest cycle.
        required_samples: usize,
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
    fn g0_periodic_surface_spectrum_is_order_invariant_and_resolved() {
        let first = PeriodicSurfaceHarmonic {
            cycles_per_track: 1,
            cosine_amplitude_m: 2.0e-6,
            sine_amplitude_m: -1.0e-6,
        };
        let third = PeriodicSurfaceHarmonic {
            cycles_per_track: 3,
            cosine_amplitude_m: 0.5e-6,
            sine_amplitude_m: 0.25e-6,
        };
        let ordered = PeriodicHarmonicSurface::new(
            "disc/material-track",
            "measured-fourier-fit",
            InputAuthority::CallerDeclared,
            0.04,
            32,
            vec![first, third],
        )
        .expect("eight samples per shortest harmonic are admitted");
        let permuted = PeriodicHarmonicSurface::new(
            "disc/material-track",
            "measured-fourier-fit",
            InputAuthority::CallerDeclared,
            0.04,
            32,
            vec![third, first],
        )
        .expect("component order is not physical");
        assert_eq!(ordered, permuted);
        assert_eq!(ordered.harmonics(), &[first, third]);

        let trace = ordered.realize().expect("bounded explicit realization");
        assert_eq!(trace.boundary(), SurfaceTraceBoundary::Periodic);
        assert_eq!(trace.heights_m().len(), 32);
        assert_eq!(trace.sample_spacing_m().to_bits(), 0.00125_f64.to_bits());
        assert!((trace.heights_m()[0] - 2.5e-6).abs() < 1.0e-20);
    }

    #[test]
    fn g0_self_affine_profile_replays_exact_power_and_requested_rms() {
        let spectrum = SelfAffinePeriodicProfileSpectrum::new(4.5e-10, 0.8, 3, 96, 0x5eed)
            .expect("bounded self-affine profile");
        let first = spectrum
            .realize_harmonics()
            .expect("deterministic Fourier realization");
        let replay = spectrum.realize_harmonics().expect("same seed must replay");
        assert_eq!(first, replay);
        assert_eq!(first.first().expect("first cycle").cycles_per_track, 3);
        assert_eq!(first.last().expect("last cycle").cycles_per_track, 96);

        let mean_square_m2 = 0.5
            * first
                .iter()
                .map(|harmonic| {
                    harmonic.cosine_amplitude_m.mul_add(
                        harmonic.cosine_amplitude_m,
                        harmonic.sine_amplitude_m * harmonic.sine_amplitude_m,
                    )
                })
                .sum::<f64>();
        let observed_rms_m = det::sqrt(mean_square_m2);
        assert!(
            (observed_rms_m / spectrum.rms_height_m() - 1.0).abs() < 2.0e-14,
            "requested RMS was not retained: requested={:.17e}, observed={observed_rms_m:.17e}",
            spectrum.rms_height_m()
        );

        let power = |harmonic: &PeriodicSurfaceHarmonic| {
            harmonic.cosine_amplitude_m.mul_add(
                harmonic.cosine_amplitude_m,
                harmonic.sine_amplitude_m * harmonic.sine_amplitude_m,
            )
        };
        let low = &first[2];
        let high = &first[31];
        let expected_ratio = det::exp(
            -(2.0 * spectrum.hurst_exponent() + 1.0)
                * det::ln(f64::from(high.cycles_per_track) / f64::from(low.cycles_per_track)),
        );
        assert!(
            (power(high) / power(low) / expected_ratio - 1.0).abs() < 2.0e-14,
            "realized Fourier power does not follow the declared self-affine PSD"
        );

        let other_seed = SelfAffinePeriodicProfileSpectrum::new(4.5e-10, 0.8, 3, 96, 0x5eee)
            .expect("second bounded realization")
            .realize_harmonics()
            .expect("second deterministic realization");
        assert_ne!(first, other_seed, "phase seed must affect geometry");
        for (a, b) in first.iter().zip(&other_seed) {
            assert!(
                (power(a) / power(b) - 1.0).abs() < 2.0e-14,
                "phase seed must not alter the declared PSD"
            );
        }
    }

    #[test]
    fn g0_self_affine_profile_refuses_invalid_statistics_and_band() {
        for hurst in [f64::NAN, 0.0, 1.0] {
            assert!(matches!(
                SelfAffinePeriodicProfileSpectrum::new(1.0e-9, hurst, 1, 8, 0),
                Err(SurfaceExcitationError::InvalidInput {
                    field: "self_affine_hurst_exponent"
                })
            ));
        }
        assert!(matches!(
            SelfAffinePeriodicProfileSpectrum::new(1.0e-9, 0.8, 0, 8, 0),
            Err(SurfaceExcitationError::InvalidInput {
                field: "self_affine_cycle_band"
            })
        ));
        assert!(matches!(
            SelfAffinePeriodicProfileSpectrum::new(1.0e-9, 0.8, 9, 8, 0),
            Err(SurfaceExcitationError::InvalidInput {
                field: "self_affine_cycle_band"
            })
        ));
        let first_unresolvable_cycle =
            u32::try_from(MAX_TRACE_SAMPLES / MIN_SAMPLES_PER_SHORTEST_HARMONIC + 1)
                .expect("trace ceiling fits u32");
        assert!(matches!(
            SelfAffinePeriodicProfileSpectrum::new(
                1.0e-9,
                0.8,
                first_unresolvable_cycle,
                first_unresolvable_cycle,
                0,
            ),
            Err(SurfaceExcitationError::SurfaceSpectrumUnderresolved {
                sample_count: MAX_TRACE_SAMPLES,
                required_samples,
            }) if required_samples == MAX_TRACE_SAMPLES + MIN_SAMPLES_PER_SHORTEST_HARMONIC
        ));
    }

    #[test]
    fn g0_periodic_surface_spectrum_refuses_aliasing_and_duplicate_cycles() {
        let harmonic = PeriodicSurfaceHarmonic {
            cycles_per_track: 5,
            cosine_amplitude_m: 1.0e-6,
            sine_amplitude_m: 0.0,
        };
        assert!(matches!(
            PeriodicHarmonicSurface::new(
                "track",
                "source",
                InputAuthority::Estimated,
                0.1,
                32,
                vec![harmonic],
            ),
            Err(SurfaceExcitationError::SurfaceSpectrumUnderresolved {
                sample_count: 32,
                required_samples: 40,
            })
        ));
        assert!(matches!(
            PeriodicHarmonicSurface::new(
                "track",
                "source",
                InputAuthority::Estimated,
                0.1,
                40,
                vec![harmonic, harmonic],
            ),
            Err(SurfaceExcitationError::DuplicateHarmonicCycle)
        ));
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
    fn g0_point_surface_pair_uses_the_same_geometry_and_material_derivative_signs() {
        let rising = trace(
            "rising-periodic",
            0.001,
            vec![0.0, 1.0e-6, 2.0e-6, 3.0e-6],
            SurfaceTraceBoundary::Periodic,
        );
        let falling = trace(
            "falling-periodic",
            0.001,
            vec![0.0, -2.0e-6, -4.0e-6, -6.0e-6],
            SurfaceTraceBoundary::Periodic,
        );
        let interface = interface();
        let receipt = evaluate_point_surface_pair(
            &interface,
            SurfaceTraceMotion {
                trace: &rising,
                path_coordinate_m: 0.0005,
                path_speed_m_per_s: 0.25,
            },
            SurfaceTraceMotion {
                trace: &falling,
                path_coordinate_m: 0.0005,
                path_speed_m_per_s: -0.5,
            },
        )
        .expect("piecewise-linear point geometry is admitted");
        assert!((receipt.combined_effective_height_m + 0.5e-6).abs() < 1.0e-21);
        assert!((receipt.filtered_surface_slopes[0] - 1.0e-3).abs() < 1.0e-15);
        assert!((receipt.filtered_surface_slopes[1] + 2.0e-3).abs() < 1.0e-15);
        assert!((receipt.combined_effective_height_rate_m_per_s - 1.25e-3).abs() < 1.0e-15);
        assert_eq!(receipt.projected_half_width_m.to_bits(), 0.0_f64.to_bits());
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

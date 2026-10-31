//! Geometric-optics propagation through smooth, isotropic gradient-index media.
//!
//! The arc-length equations are `r' = t`,
//! `t' = (grad(n) - t (t · grad(n))) / n`, obtained from
//! `d(n t)/ds = grad(n)`. The phase index bends the ray; the separately supplied
//! group index integrates narrowband delay. A transverse axis follows
//! `e' = -(e · t') t`, the parallel-transport law for polarization. Power
//! attenuation is `exp(-integral sigma_a ds)` along the curved physical path.
//!
//! The allocation-free Dormand–Prince 5(4) integrator controls local estimates
//! for position, direction, transverse axis and both path integrals. These are
//! numerical estimates, not interval certificates or global error bounds.
//! Callers must resolve field features with `max_step_m`, provide a smooth
//! field and terminate segments at material discontinuities. [`event`] adds
//! bounded scene-boundary advancement; [`bundle`] adds analytic angular
//! variations and endpoint-plane focusing. Neither supplies scalar-camera/BDPT
//! MIS, scattering, birefringence, diffraction or spin-Hall corrections.
//!
//! References: Moore, JOSA 65, 451–455 (1975),
//! <https://doi.org/10.1364/JOSA.65.000451>; Dormand and Prince, JCAM 6, 19–26
//! (1980), <https://doi.org/10.1016/0771-050X(80)90013-3>; Bliokh, J. Opt. A
//! 11, 094009 (2009), <https://doi.org/10.1088/1464-4258/11/9/094009>.

use fs_geom::{Aabb, Point3, Vec3};
use fs_math::det;

use crate::charts::Ray;
use crate::dielectric::{DielectricError, DielectricGlass, DielectricIor};
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN};

pub mod bundle;
pub mod event;

/// Version of the deterministic gradient-index propagation arithmetic.
pub const GRADIENT_INDEX_SEMANTICS_VERSION: u32 = 2;
/// Hard ceiling on attempted adaptive steps, including rejected steps.
pub const MAX_GRADIENT_INDEX_STEPS: u32 = 1_000_000;

/// Admission or execution refusal. Propagation never publishes a partial result.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GradientIndexError {
    /// A phase index, gradient, absorption or group index is inadmissible.
    InvalidSample,
    /// A position is outside the caller's explicitly supported field domain.
    OutsideFieldDomain,
    /// The ray has a nonfinite position or a zero/nonfinite direction.
    InvalidRay,
    /// The transverse reference is zero, nonfinite or parallel to the ray.
    InvalidTransverseAxis,
    /// The requested distance or step/tolerance/work envelope is invalid.
    InvalidOptions,
    /// A material's spectral owner refused the query.
    Dielectric(DielectricError),
    /// A group index was present initially but absent later on the path.
    GroupIndexUnavailable,
    /// A finite result cannot be represented or a step cannot advance.
    NumericalRange,
    /// The local estimate exceeds tolerance at the minimum admitted step.
    AccuracyNotReached,
    /// The maximum attempted-step count was consumed before completing the path.
    StepBudgetExceeded,
    /// The caller's bounded-work control requested cancellation.
    Cancelled,
    /// A geometric clearance/proximity query was malformed or contradictory.
    InvalidBoundaryQuery,
    /// A nearby surface lacked a resolvable transverse contact or safe progress.
    BoundaryResolutionLimit,
    /// The previous contact has no admitted local departure neighborhood.
    UnsupportedBoundaryDeparture,
    /// The explicit per-path geometry-query budget was exhausted.
    BoundaryQueryBudgetExceeded,
    /// A trial stage left its geometric free ball or admitted departure side.
    /// Event propagation retries with a shorter step; public fixed-distance
    /// propagation does not produce this refusal.
    BoundaryClearanceExceeded,
}

impl From<DielectricError> for GradientIndexError {
    fn from(value: DielectricError) -> Self {
        Self::Dielectric(value)
    }
}

impl core::fmt::Display for GradientIndexError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "gradient-index propagation: {self:?}")
    }
}

impl core::error::Error for GradientIndexError {}

/// Optical properties and their analytic spatial derivative at one wavelength.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GradientIndexSample {
    phase_index: f64,
    phase_gradient_per_m: Vec3,
    absorption_per_m: f64,
    group_index: Option<f64>,
}

impl GradientIndexSample {
    /// Admit an absolute phase index in `[1, 3]`, its finite gradient in inverse
    /// metres, a nonnegative power-absorption coefficient in inverse metres,
    /// and an optional positive finite group index. The group index must come
    /// from spectral differentiation; it is never inferred from phase index.
    pub fn try_new(
        phase_index: f64,
        phase_gradient_per_m: Vec3,
        absorption_per_m: f64,
        group_index: Option<f64>,
    ) -> Result<Self, GradientIndexError> {
        if !phase_index.is_finite()
            || !(1.0..=3.0).contains(&phase_index)
            || !finite_vec(phase_gradient_per_m)
            || !absorption_per_m.is_finite()
            || absorption_per_m < 0.0
            || group_index.is_some_and(|n| !n.is_finite() || n <= 0.0)
        {
            return Err(GradientIndexError::InvalidSample);
        }
        Ok(Self {
            phase_index,
            phase_gradient_per_m,
            absorption_per_m,
            group_index,
        })
    }

    /// Absolute phase index at the queried vacuum wavelength.
    #[must_use]
    pub const fn phase_index(self) -> f64 {
        self.phase_index
    }

    /// Spatial derivative of phase index, in inverse metres.
    #[must_use]
    pub const fn phase_gradient_per_m(self) -> Vec3 {
        self.phase_gradient_per_m
    }

    /// Beer–Lambert power-absorption coefficient, in inverse metres.
    #[must_use]
    pub const fn absorption_per_m(self) -> f64 {
        self.absorption_per_m
    }

    /// Narrowband group index, if spectral differentiation is supported here.
    #[must_use]
    pub const fn group_index(self) -> Option<f64> {
        self.group_index
    }
}

/// A spatial optical law pinned to one immutable physical-state snapshot.
///
/// Positions are metres and wavelengths are vacuum nanometres. The law must
/// return consistent phase values and derivatives, remain smooth over each
/// propagated segment and refuse queries outside its physical/source domain.
/// State-dependent adapters borrow or own the accepted state; sampling must
/// not silently refresh it while a path is being integrated.
pub trait GradientIndexField {
    /// Evaluate the local spectral optical law and analytic spatial gradient.
    fn sample(
        &self,
        position: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, GradientIndexError>;

    /// Whether the optical law is spatially constant over all finite positions.
    /// The default keeps bounded or varying fields on the adaptive path. A
    /// `true` declaration permits exact straight propagation from one sample;
    /// it must never be inferred from a gradient that vanishes at one point.
    fn is_spatially_uniform(&self) -> bool {
        false
    }
}

/// Borrowed native spectral glass with no spatial refractive-index gradient.
#[derive(Clone, Copy, Debug)]
pub struct HomogeneousIndexField<'a> {
    glass: &'a DielectricGlass,
}

impl<'a> HomogeneousIndexField<'a> {
    /// Reuse the exact native dispersion and absorption laws without copying tables.
    #[must_use]
    pub const fn new(glass: &'a DielectricGlass) -> Self {
        Self { glass }
    }
}

impl GradientIndexField for HomogeneousIndexField<'_> {
    fn sample(
        &self,
        position: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        if !finite_point(position) {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        material_sample(self.glass, wavelength_nm, 0.0, Vec3::new(0.0, 0.0, 0.0))
    }

    fn is_spatially_uniform(&self) -> bool {
        true
    }
}

/// A bounded affine spatial index on top of a native spectral dispersion law:
/// `n(r, lambda) = n_reference(lambda) + gradient · (r - reference_point)`.
///
/// The spatial offset is wavelength-independent. Therefore it contributes the
/// same offset to group index, while absorption retains the borrowed spectral
/// material law. This is an explicit model, not an inferred temperature law.
#[derive(Clone, Copy, Debug)]
pub struct AffineIndexField<'a> {
    glass: &'a DielectricGlass,
    reference_point: Point3,
    gradient_per_m: Vec3,
    bounds: Aabb,
}

impl<'a> AffineIndexField<'a> {
    /// Admit finite bounds and keep the complete visible phase-index law in
    /// `[1, 3]` throughout the box. Native analytic laws are monotone on their
    /// admitted visible domain; every tabulated knot is checked separately.
    pub fn try_new(
        glass: &'a DielectricGlass,
        reference_point: Point3,
        gradient_per_m: Vec3,
        bounds: Aabb,
    ) -> Result<Self, GradientIndexError> {
        if !finite_point(reference_point)
            || !finite_vec(gradient_per_m)
            || !finite_point(bounds.min)
            || !finite_point(bounds.max)
            || bounds.min.x >= bounds.max.x
            || bounds.min.y >= bounds.max.y
            || bounds.min.z >= bounds.max.z
        {
            return Err(GradientIndexError::InvalidSample);
        }
        let candidate = Self {
            glass,
            reference_point,
            gradient_per_m,
            bounds,
        };
        for x in [bounds.min.x, bounds.max.x] {
            for y in [bounds.min.y, bounds.max.y] {
                for z in [bounds.min.z, bounds.max.z] {
                    let p = Point3::new(x, y, z);
                    candidate.sample(p, LAMBDA_MIN)?;
                    candidate.sample(p, LAMBDA_MAX)?;
                    if let DielectricIor::Tabulated(table) = glass.ior() {
                        for knot in table.samples() {
                            candidate.sample(p, knot.wavelength_nm())?;
                        }
                    }
                }
            }
        }
        Ok(candidate)
    }
}

impl GradientIndexField for AffineIndexField<'_> {
    fn sample(
        &self,
        position: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        if !finite_point(position) || !self.bounds.contains(position) {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        let offset = self
            .gradient_per_m
            .dot(position.delta_from(self.reference_point));
        material_sample(self.glass, wavelength_nm, offset, self.gradient_per_m)
    }
}

fn material_sample(
    glass: &DielectricGlass,
    wavelength_nm: f64,
    offset: f64,
    gradient: Vec3,
) -> Result<GradientIndexSample, GradientIndexError> {
    let n = glass.ior().eval(wavelength_nm)? + offset;
    let group = match glass.ior().group_index(wavelength_nm) {
        Ok(n) if (n + offset).is_finite() && n + offset > 0.0 => Some(n + offset),
        // A steep supplied phase interpolant may imply a nonpositive group
        // velocity model. Phase refraction remains usable, but this owner
        // does not claim anomalous-dispersion pulse delays from that table.
        Ok(_) => None,
        Err(DielectricError::UndefinedGroupIndex) => None,
        Err(error) => return Err(error.into()),
    };
    GradientIndexSample::try_new(
        n,
        gradient,
        glass.absorption_ref().extinction_per_m(wavelength_nm)?,
        group,
    )
}

/// Local error controls and explicit metric step/work bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GradientIndexOptions {
    /// Absolute local position-estimate tolerance, metres.
    pub position_tolerance_m: f64,
    /// Absolute local unit-direction and transverse-axis tolerance.
    pub direction_tolerance: f64,
    /// Absolute local phase/group optical-path tolerance, metres.
    pub optical_path_tolerance_m: f64,
    /// Absolute local dimensionless absorption-depth tolerance.
    pub optical_depth_tolerance: f64,
    /// Relative tolerance on each local increment (never on global coordinates).
    pub relative_tolerance: f64,
    /// Initial trial length in metres, clamped to the remaining distance.
    pub initial_step_m: f64,
    /// Smallest retry length; a final shortened endpoint step may be smaller.
    pub min_step_m: f64,
    /// Largest trial length; caller must use this to resolve narrow field features.
    pub max_step_m: f64,
    /// Maximum attempted steps, including rejected trials.
    pub max_steps: u32,
}

impl Default for GradientIndexOptions {
    fn default() -> Self {
        Self {
            position_tolerance_m: 1.0e-9,
            direction_tolerance: 1.0e-10,
            optical_path_tolerance_m: 1.0e-9,
            optical_depth_tolerance: 1.0e-10,
            relative_tolerance: 1.0e-9,
            initial_step_m: 0.01,
            min_step_m: 1.0e-10,
            max_step_m: 0.1,
            max_steps: 16_384,
        }
    }
}

impl GradientIndexOptions {
    fn validate(self) -> Result<(), GradientIndexError> {
        if [
            self.position_tolerance_m,
            self.direction_tolerance,
            self.optical_path_tolerance_m,
            self.optical_depth_tolerance,
            self.initial_step_m,
            self.min_step_m,
            self.max_step_m,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v <= 0.0)
            || !self.relative_tolerance.is_finite()
            || self.relative_tolerance < 0.0
            || self.min_step_m > self.initial_step_m
            || self.initial_step_m > self.max_step_m
            || self.max_steps == 0
            || self.max_steps > MAX_GRADIENT_INDEX_STEPS
        {
            return Err(GradientIndexError::InvalidOptions);
        }
        Ok(())
    }
}

/// Completed curved propagation, including integrals along physical arc length.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GradientIndexPropagation {
    /// Endpoint position and normalized tangent; chart-ray parameter is metres.
    pub ray: Ray,
    /// Parallel-transported normalized transverse reference axis.
    pub transverse_axis: Vec3,
    /// Exact field query at the initial ray origin.
    pub initial_sample: GradientIndexSample,
    /// Local optical law at the returned endpoint position. Uniform fields
    /// reuse their constant initial sample without another spectral evaluation.
    pub final_sample: GradientIndexSample,
    /// Completed geometric arc length, metres.
    pub distance_m: f64,
    /// Integral of phase index over arc length, metres.
    pub optical_path_length_m: f64,
    /// Integral of group index over arc length, metres; divide by c for delay.
    /// Absent when the initial field sample has no supported group index.
    pub group_optical_path_length_m: Option<f64>,
    /// Integral of power-absorption coefficient over arc length.
    pub optical_depth: f64,
    /// Number of accepted integration steps (one for a nonzero uniform segment).
    pub accepted_steps: u32,
    /// Number of rejected adaptive steps.
    pub rejected_steps: u32,
    /// Total field queries, including rejected trials.
    pub field_evaluations: u32,
    /// Largest normalized local error estimate on an accepted step, at most 1.
    pub max_normalized_error_estimate: f64,
}

impl GradientIndexPropagation {
    /// Beam power transmission. This does not include interface Fresnel loss.
    #[must_use]
    pub fn power_transmittance(self) -> f64 {
        det::exp(-self.optical_depth)
    }

    /// Absorbed power fraction, retaining optically thin losses without `1-exp`.
    #[must_use]
    pub fn absorbed_power_fraction(self) -> f64 {
        -det::expm1(-self.optical_depth)
    }
}

/// Propagate a ray and transverse frame through one smooth optical field segment.
pub fn propagate<F: GradientIndexField + ?Sized>(
    field: &F,
    ray: Ray,
    transverse_axis: Vec3,
    wavelength_nm: f64,
    distance_m: f64,
    options: GradientIndexOptions,
) -> Result<GradientIndexPropagation, GradientIndexError> {
    propagate_with_control(
        field,
        ray,
        transverse_axis,
        wavelength_nm,
        distance_m,
        options,
        || true,
    )
}

/// Controlled propagation; `control` is checked before initial sampling and
/// every attempted step, including consecutive rejections. Returning `false`
/// refuses the operation with `Cancelled`. A caller can poll its `Cx` here.
/// No callback is invoked after the returned complete result is constructed.
pub fn propagate_with_control<F: GradientIndexField + ?Sized>(
    field: &F,
    ray: Ray,
    transverse_axis: Vec3,
    wavelength_nm: f64,
    distance_m: f64,
    options: GradientIndexOptions,
    mut control: impl FnMut() -> bool,
) -> Result<GradientIndexPropagation, GradientIndexError> {
    options.validate()?;
    if !distance_m.is_finite() || distance_m < 0.0 {
        return Err(GradientIndexError::InvalidOptions);
    }
    if !wavelength_nm.is_finite() || !(LAMBDA_MIN..=LAMBDA_MAX).contains(&wavelength_nm) {
        return Err(DielectricError::InvalidWavelength.into());
    }
    if !finite_point(ray.origin) {
        return Err(GradientIndexError::InvalidRay);
    }
    let direction = unit(ray.dir).ok_or(GradientIndexError::InvalidRay)?;
    let axis =
        transverse(direction, transverse_axis).ok_or(GradientIndexError::InvalidTransverseAxis)?;
    if !control() {
        return Err(GradientIndexError::Cancelled);
    }
    let initial_sample = field.sample(ray.origin, wavelength_nm)?;
    if field.is_spatially_uniform() {
        if initial_sample.phase_gradient_per_m != Vec3::new(0.0, 0.0, 0.0) {
            return Err(GradientIndexError::InvalidSample);
        }
        let origin = ray.origin.offset(direction.scale(distance_m));
        let optical_path_length_m = initial_sample.phase_index * distance_m;
        let group_optical_path_length_m = initial_sample.group_index.map(|n| n * distance_m);
        let optical_depth = initial_sample.absorption_per_m * distance_m;
        if !finite_point(origin)
            || !optical_path_length_m.is_finite()
            || group_optical_path_length_m.is_some_and(|v| !v.is_finite())
            || !optical_depth.is_finite()
        {
            return Err(GradientIndexError::NumericalRange);
        }
        return Ok(GradientIndexPropagation {
            ray: Ray {
                origin,
                dir: direction,
            },
            transverse_axis: axis,
            initial_sample,
            final_sample: initial_sample,
            distance_m,
            optical_path_length_m,
            group_optical_path_length_m,
            optical_depth,
            accepted_steps: u32::from(distance_m > 0.0),
            rejected_steps: 0,
            field_evaluations: 1,
            max_normalized_error_estimate: 0.0,
        });
    }
    let track_group = initial_sample.group_index.is_some();
    // Integrate displacement from the ray origin: local tolerances must not
    // inflate when an otherwise identical scene is translated far from zero.
    let mut state = [0.0; 12];
    put_vec(&mut state, 3, direction);
    put_vec(&mut state, 6, axis);
    let mut final_sample = initial_sample;
    let mut derivative = rhs(&state, initial_sample, track_group)?;
    let mut distance = 0.0;
    let mut step = options.initial_step_m.min(distance_m);
    let mut accepted = 0;
    let mut rejected = 0;
    let mut evaluations = 1;
    let mut max_error = 0.0_f64;

    while distance < distance_m {
        if !control() {
            return Err(GradientIndexError::Cancelled);
        }
        if accepted + rejected >= options.max_steps {
            return Err(GradientIndexError::StepBudgetExceeded);
        }
        let remaining = distance_m - distance;
        step = step.min(remaining);
        if !step.is_finite() || step <= 0.0 || distance + step == distance {
            return Err(GradientIndexError::NumericalRange);
        }
        let trial = trial_step(
            field,
            ray.origin,
            wavelength_nm,
            &state,
            derivative,
            step,
            track_group,
        )?;
        evaluations += 6;
        let error = normalized_error(&state, &trial, step, options)?;
        let factor = if error == 0.0 {
            5.0
        } else {
            (0.9 * det::pow(error, -0.2)).clamp(0.2, 5.0)
        };
        if error <= 1.0 {
            state = trial.state;
            final_sample = trial.sample;
            derivative = rhs(&state, final_sample, track_group)?;
            distance = if step == remaining {
                distance_m
            } else {
                distance + step
            };
            accepted += 1;
            max_error = max_error.max(error);
            step = (step * factor)
                .min(options.max_step_m)
                .max(options.min_step_m);
        } else {
            rejected += 1;
            if step <= options.min_step_m {
                return Err(GradientIndexError::AccuracyNotReached);
            }
            step = (step * factor.min(0.9)).max(options.min_step_m);
        }
    }
    let origin = ray.origin.offset(get_vec(&state, 0));
    if !finite_point(origin) {
        return Err(GradientIndexError::NumericalRange);
    }
    Ok(GradientIndexPropagation {
        ray: Ray {
            origin,
            dir: get_vec(&state, 3),
        },
        transverse_axis: get_vec(&state, 6),
        initial_sample,
        final_sample,
        distance_m,
        optical_path_length_m: state[9],
        group_optical_path_length_m: track_group.then_some(state[10]),
        optical_depth: state[11],
        accepted_steps: accepted,
        rejected_steps: rejected,
        field_evaluations: evaluations,
        max_normalized_error_estimate: max_error,
    })
}

type State = [f64; 12];

struct Trial {
    state: State,
    error: State,
    projection_error: f64,
    sample: GradientIndexSample,
}

// Dormand–Prince 5(4), with the first derivative reused across rejected trials.
// The endpoint field query is reused after accepted tangent/frame projection,
// so every trial costs six new field evaluations rather than seven.
const A: [[f64; 6]; 6] = [
    [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
    [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
    [
        19372.0 / 6561.0,
        -25360.0 / 2187.0,
        64448.0 / 6561.0,
        -212.0 / 729.0,
        0.0,
        0.0,
    ],
    [
        9017.0 / 3168.0,
        -355.0 / 33.0,
        46732.0 / 5247.0,
        49.0 / 176.0,
        -5103.0 / 18656.0,
        0.0,
    ],
    [
        35.0 / 384.0,
        0.0,
        500.0 / 1113.0,
        125.0 / 192.0,
        -2187.0 / 6784.0,
        11.0 / 84.0,
    ],
];
const ERROR: [f64; 7] = [
    35.0 / 384.0 - 5179.0 / 57600.0,
    0.0,
    500.0 / 1113.0 - 7571.0 / 16695.0,
    125.0 / 192.0 - 393.0 / 640.0,
    -2187.0 / 6784.0 + 92097.0 / 339200.0,
    11.0 / 84.0 - 187.0 / 2100.0,
    -1.0 / 40.0,
];

fn trial_step<F: GradientIndexField + ?Sized>(
    field: &F,
    origin: Point3,
    wavelength_nm: f64,
    state: &State,
    first_derivative: State,
    step: f64,
    track_group: bool,
) -> Result<Trial, GradientIndexError> {
    trial_step_guarded(
        field,
        origin,
        wavelength_nm,
        state,
        first_derivative,
        step,
        track_group,
        |_, _| true,
        &mut 0,
    )
}

#[allow(clippy::too_many_arguments)]
fn trial_step_guarded<F: GradientIndexField + ?Sized>(
    field: &F,
    origin: Point3,
    wavelength_nm: f64,
    state: &State,
    first_derivative: State,
    step: f64,
    track_group: bool,
    mut admits_stage: impl FnMut(Point3, Vec3) -> bool,
    field_evaluations: &mut u32,
) -> Result<Trial, GradientIndexError> {
    let mut k = [[0.0; 12]; 7];
    k[0] = first_derivative;
    let mut stage = *state;
    let mut last_sample = None;
    for (i, row) in A.iter().enumerate() {
        for j in 0..12 {
            let mut increment = 0.0;
            for a in 0..=i {
                increment += row[a] * k[a][j];
            }
            stage[j] = state[j] + step * increment;
        }
        if stage.iter().any(|x| !x.is_finite()) {
            return Err(GradientIndexError::NumericalRange);
        }
        let point = origin.offset(get_vec(&stage, 0));
        if !finite_point(point) {
            return Err(GradientIndexError::NumericalRange);
        }
        if !admits_stage(point, get_vec(&stage, 3)) {
            return Err(GradientIndexError::BoundaryClearanceExceeded);
        }
        *field_evaluations += 1;
        let sample = field.sample(point, wavelength_nm)?;
        k[i + 1] = rhs(&stage, sample, track_group)?;
        last_sample = Some(sample);
    }
    let mut error = [0.0; 12];
    for j in 0..12 {
        for i in 0..7 {
            error[j] += ERROR[i] * k[i][j];
        }
        error[j] *= step;
    }
    let raw_direction = get_vec(&stage, 3);
    let raw_axis = get_vec(&stage, 6);
    let direction = unit(raw_direction).ok_or(GradientIndexError::NumericalRange)?;
    let axis = transverse(direction, raw_axis).ok_or(GradientIndexError::NumericalRange)?;
    // Projection restores the analytic unit/orthogonality invariants; include
    // its correction in the local estimate so it cannot conceal a large step.
    let projection_error = norm(sub(direction, raw_direction)).max(norm(sub(axis, raw_axis)));
    put_vec(&mut stage, 3, direction);
    put_vec(&mut stage, 6, axis);
    if stage[9] < state[9] || stage[10] < state[10] || stage[11] < state[11] {
        return Err(GradientIndexError::NumericalRange);
    }
    Ok(Trial {
        state: stage,
        error,
        projection_error,
        sample: last_sample.ok_or(GradientIndexError::NumericalRange)?,
    })
}

fn rhs(
    state: &State,
    sample: GradientIndexSample,
    track_group: bool,
) -> Result<State, GradientIndexError> {
    let direction = unit(get_vec(state, 3)).ok_or(GradientIndexError::NumericalRange)?;
    let gradient = sample.phase_gradient_per_m.scale(1.0 / sample.phase_index);
    let curvature = sub(gradient, direction.scale(direction.dot(gradient)));
    let axis_derivative = direction.scale(-get_vec(state, 6).dot(curvature));
    let mut result = [0.0; 12];
    put_vec(&mut result, 0, direction);
    put_vec(&mut result, 3, curvature);
    put_vec(&mut result, 6, axis_derivative);
    result[9] = sample.phase_index;
    result[10] = if track_group {
        sample
            .group_index
            .ok_or(GradientIndexError::GroupIndexUnavailable)?
    } else {
        0.0
    };
    result[11] = sample.absorption_per_m;
    if result.iter().any(|x| !x.is_finite()) {
        return Err(GradientIndexError::NumericalRange);
    }
    Ok(result)
}

fn normalized_error(
    state: &State,
    trial: &Trial,
    step: f64,
    options: GradientIndexOptions,
) -> Result<f64, GradientIndexError> {
    let angular_scale = options.direction_tolerance + options.relative_tolerance;
    let mut result = trial.projection_error / angular_scale;
    for i in 0..12 {
        let scale = match i {
            0..=2 => options.position_tolerance_m + options.relative_tolerance * step,
            3..=8 => angular_scale,
            9..=10 => {
                options.optical_path_tolerance_m
                    + options.relative_tolerance * (trial.state[i] - state[i]).abs()
            }
            _ => {
                options.optical_depth_tolerance
                    + options.relative_tolerance * (trial.state[i] - state[i]).abs()
            }
        };
        let error = trial.error[i].abs() / scale;
        if !scale.is_finite() || scale <= 0.0 || !error.is_finite() {
            return Err(GradientIndexError::NumericalRange);
        }
        result = result.max(error);
    }
    if !result.is_finite() {
        return Err(GradientIndexError::NumericalRange);
    }
    Ok(result)
}

fn finite_point(p: Point3) -> bool {
    p.x.is_finite() && p.y.is_finite() && p.z.is_finite()
}
fn finite_vec(v: Vec3) -> bool {
    v.x.is_finite() && v.y.is_finite() && v.z.is_finite()
}
fn get_vec(state: &State, offset: usize) -> Vec3 {
    Vec3::new(state[offset], state[offset + 1], state[offset + 2])
}
fn put_vec(state: &mut State, offset: usize, v: Vec3) {
    state[offset] = v.x;
    state[offset + 1] = v.y;
    state[offset + 2] = v.z;
}
fn sub(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}
fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}
fn norm(v: Vec3) -> f64 {
    det::sqrt(v.dot(v))
}
fn unit(v: Vec3) -> Option<Vec3> {
    if !finite_vec(v) {
        return None;
    }
    let scale = v.x.abs().max(v.y.abs()).max(v.z.abs());
    if scale == 0.0 {
        return None;
    }
    let scaled = Vec3::new(v.x / scale, v.y / scale, v.z / scale);
    Some(scaled.scale(1.0 / norm(scaled)))
}
fn transverse(direction: Vec3, reference: Vec3) -> Option<Vec3> {
    let reference = unit(reference)?;
    let vertical = unit(cross(direction, reference))?;
    unit(cross(vertical, direction))
}

#[cfg(test)]
mod tests;

//! Shear-sensitive Findley fatigue-limit usage for a resolved closed stress block.
//!
//! On each declared fixed plane, the shear amplitude is the radius of the
//! smallest circle enclosing the two-dimensional shear-traction path. The
//! selected plane maximizes `tau_amplitude + k * maximum_normal_stress` over the
//! ENTIRE supplied block; shear and normal extrema need not coincide. Usage is
//! `max(0, selected_parameter) / f`, including the named zero clipping under
//! sufficiently compressive normal stress. k and f must be supplied as a
//! transferable Findley calibration, not inferred from static strength or an
//! axial fatigue coefficient. See the Findley and critical-plane definitions:
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.fatigue/fatigue_ug_sme.4.22.html>
//! <https://doc.comsol.com/6.3/doc/com.comsol.help.fatigue/fatigue_ug_sme.4.18.html>
//!
//! This is a fatigue-LIMIT comparison for a caller-declared representative
//! periodic, elastic stress block. It does not count its internal reversals,
//! infer repetitions, predict finite life, accumulate Miner damage, or establish
//! cyclic stabilization or material applicability. Every physical extremum must
//! be resolved by an observation. Closure is checked on the FULL tensor; a
//! constant tensor history is refused rather than treated as a loading cycle.
//! All stresses (Pa) and plane normals share one fixed frame, in tensor-Voigt
//! order `[xx, yy, zz, xy, yz, zx]`.
//!
//! Tensor history is stored once; planes have precomputed traction projections.
//! Finalization reuses one O(samples) shear-path buffer, so total storage is
//! O(samples + planes). A normalized, translated, deterministic incremental
//! enclosing-circle solver has an explicit work cap and verifies enclosure.
//! This floating-point calculation and the finite plane search are numerical
//! model estimates, not exact-geometry certificates or a continuum optimum.

pub mod source;

use crate::fatigue::{FatigueError, FatigueScope, NormalStressDirection};
use fs_math::det;
use std::sync::Arc;

/// Already calibrated material parameters, with their declared applicability.
#[derive(Debug, Clone, PartialEq)]
pub struct FindleyCalibration {
    scope: FatigueScope,
    normal_stress_sensitivity: f64,
    limit_factor_pa: f64,
    source: Option<Arc<source::FindleyCalibrationSource>>,
}
impl FindleyCalibration {
    pub fn new(
        scope: FatigueScope,
        normal_stress_sensitivity: f64,
        limit_factor_pa: f64,
    ) -> Result<Self, FindleyError> {
        scope.validate()?;
        if !normal_stress_sensitivity.is_finite() || normal_stress_sensitivity < 0.0 {
            return Err(FindleyError::Invalid(
                "Findley k must be finite and nonnegative",
            ));
        }
        if !limit_factor_pa.is_finite() || limit_factor_pa <= 0.0 {
            return Err(FindleyError::Invalid(
                "Findley f must be a finite positive pressure",
            ));
        }
        Ok(Self {
            scope,
            normal_stress_sensitivity,
            limit_factor_pa,
            source: None,
        })
    }
    pub fn scope(&self) -> &FatigueScope {
        &self.scope
    }
    pub fn normal_stress_sensitivity(&self) -> f64 {
        self.normal_stress_sensitivity
    }
    pub fn limit_factor_pa(&self) -> f64 {
        self.limit_factor_pa
    }
    pub fn source(&self) -> Option<&source::FindleyCalibrationSource> {
        self.source.as_deref()
    }
}

/// Projection work is samples times planes, admitted before accepting history.
/// A circle work unit is a 2D point normalization, permutation, containment
/// predicate, circle construction, or final enclosure check. The incremental
/// algorithm has a cubic worst case, bounded by this explicit operation cap.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FindleyLimits {
    pub maximum_samples: u64,
    pub maximum_planes: usize,
    pub maximum_projection_work: u64,
    pub maximum_circle_work: u64,
}

/// Componentwise endpoint tolerance in the ORIGINAL stress tensor:
/// `abs(last-first) <= absolute_stress_pa + relative * max(abs(first),abs(last))`.
/// No observation is modified to meet this caller-declared tolerance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StressClosureTolerance {
    pub absolute_stress_pa: f64,
    pub relative: f64,
}

/// Fatigue-limit usage on the declared finite plane set, not cycles or damage.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FindleyCycleResult {
    pub samples: u64,
    pub evaluated_planes: usize,
    pub selected_plane: usize,
    pub normal_in_frame: [f64; 3],
    pub shear_amplitude_pa: f64,
    pub minimum_normal_stress_pa: f64,
    pub maximum_normal_stress_pa: f64,
    /// Signed value before the Findley compression clipping convention.
    pub unclipped_findley_parameter_pa: f64,
    /// `max(0, tau_amplitude + k * maximum_normal_stress)`.
    pub findley_parameter_pa: f64,
    pub usage_factor: f64,
    pub projection_work: u64,
    pub circle_work: u64,
    pub closure: StressClosureTolerance,
    pub maximum_stress_closure_defect_pa: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum FindleyError {
    Invalid(&'static str),
    Budget,
    InvalidSample {
        sample: u64,
        component: usize,
    },
    NotClosed {
        component: usize,
    },
    NotCyclic,
    /// A finite, enclosed circle could not be established numerically.
    NumericalCircle,
    Fatigue(FatigueError),
    Allocation,
    Cancelled,
    Finished,
}
impl core::fmt::Display for FindleyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Findley stress-block evaluation refused: {self:?}")
    }
}
impl std::error::Error for FindleyError {}
impl From<FatigueError> for FindleyError {
    fn from(error: FatigueError) -> Self {
        match error {
            FatigueError::Cancelled => Self::Cancelled,
            FatigueError::Allocation => Self::Allocation,
            error => Self::Fatigue(error),
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct SearchPlane {
    normal: [f64; 3],
    normal_weights: [f64; 6],
    shear_weights: [[f64; 6]; 2],
}
impl SearchPlane {
    fn new(direction: [f64; 3]) -> Result<Self, FindleyError> {
        let normal = NormalStressDirection::new(direction)?.unit_normal();
        // Cross with the least aligned cardinal axis to avoid a small norm.
        let mut axis = 0;
        for j in 1..3 {
            if normal[j].abs() < normal[axis].abs() {
                axis = j;
            }
        }
        let mut cardinal = [0.0; 3];
        cardinal[axis] = 1.0;
        let tangent = NormalStressDirection::new(cross3(normal, cardinal))?.unit_normal();
        let bitangent = NormalStressDirection::new(cross3(normal, tangent))?.unit_normal();
        Ok(Self {
            normal,
            normal_weights: traction_weights(normal, normal),
            shear_weights: [
                traction_weights(tangent, normal),
                traction_weights(bitangent, normal),
            ],
        })
    }
}
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn traction_weights(a: [f64; 3], n: [f64; 3]) -> [f64; 6] {
    [
        a[0] * n[0],
        a[1] * n[1],
        a[2] * n[2],
        a[0] * n[1] + a[1] * n[0],
        a[1] * n[2] + a[2] * n[1],
        a[2] * n[0] + a[0] * n[2],
    ]
}
fn dot6(a: [f64; 6], b: [f64; 6]) -> f64 {
    let mut result = 0.0;
    for i in 0..6 {
        result += a[i] * b[i];
    }
    result
}
fn normal_stress(weights: [f64; 6], stress: [f64; 6]) -> Result<f64, FindleyError> {
    // Most physical projections fit directly. Besides avoiding six divisions,
    // this preserves a small resolved component next to a large unrelated one.
    let direct = dot6(weights, stress);
    if direct.is_finite() {
        return Ok(direct);
    }
    // Scaling is only needed when the direct sum overflows before cancellation.
    let scale = stress.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    if scale == 0.0 {
        return Ok(0.0);
    }
    let value = dot6(weights, stress.map(|x| x / scale)) * scale;
    if !value.is_finite() {
        return Err(FatigueError::Arithmetic.into());
    }
    Ok(value)
}

/// Transactional observations and finalization. Failed pushes leave history
/// unchanged; failed finalization publishes no result and may be retried.
/// Projection/circle work is performed at finish, after bounded history capture.
#[derive(Debug, Clone)]
pub struct FindleyCycle {
    calibration: Arc<FindleyCalibration>,
    planes: Arc<Vec<SearchPlane>>,
    limits: FindleyLimits,
    closure: StressClosureTolerance,
    samples: Vec<[f64; 6]>,
    result: Option<FindleyCycleResult>,
}
impl FindleyCycle {
    pub fn new(
        calibration: Arc<FindleyCalibration>,
        normals: Vec<[f64; 3]>,
        limits: FindleyLimits,
        closure: StressClosureTolerance,
    ) -> Result<Self, FindleyError> {
        if limits.maximum_samples < 3
            || limits.maximum_planes == 0
            || limits.maximum_projection_work == 0
            || limits.maximum_circle_work == 0
        {
            return Err(FindleyError::Invalid(
                "positive work/plane limits and at least three sample slots are required",
            ));
        }
        if !closure.absolute_stress_pa.is_finite()
            || closure.absolute_stress_pa < 0.0
            || !closure.relative.is_finite()
            || closure.relative < 0.0
            || closure.relative > 1.0
        {
            return Err(FindleyError::Invalid(
                "finite nonnegative closure tolerances and relative tolerance at most one are required",
            ));
        }
        if normals.is_empty() {
            return Err(FindleyError::Invalid("declare at least one search plane"));
        }
        if normals.len() > limits.maximum_planes {
            return Err(FindleyError::Budget);
        }
        let count = u64::try_from(normals.len()).map_err(|_| FindleyError::Budget)?;
        if count > limits.maximum_projection_work {
            return Err(FindleyError::Budget);
        }
        let mut planes = Vec::new();
        planes
            .try_reserve_exact(normals.len())
            .map_err(|_| FindleyError::Allocation)?;
        for normal in normals {
            planes.push(SearchPlane::new(normal)?);
        }
        Ok(Self {
            calibration,
            planes: Arc::new(planes),
            limits,
            closure,
            samples: Vec::new(),
            result: None,
        })
    }
    pub fn calibration(&self) -> &FindleyCalibration {
        &self.calibration
    }
    pub fn sample_count(&self) -> u64 {
        self.samples.len() as u64
    }
    pub fn plane_count(&self) -> usize {
        self.planes.len()
    }
    pub fn projection_work(&self) -> u64 {
        self.sample_count() * self.planes.len() as u64
    }
    pub fn circle_work(&self) -> u64 {
        self.result.map_or(0, |r| r.circle_work)
    }
    pub fn limits(&self) -> FindleyLimits {
        self.limits
    }
    pub fn closure_tolerance(&self) -> StressClosureTolerance {
        self.closure
    }
    pub fn result(&self) -> Option<&FindleyCycleResult> {
        self.result.as_ref()
    }

    /// Empty window with the SAME calibration, search planes and limits. The
    /// caller explicitly supplies its new first sample, including any boundary.
    pub fn fresh(&self) -> Result<Self, FindleyError> {
        Ok(Self {
            calibration: self.calibration.clone(),
            planes: self.planes.clone(),
            limits: self.limits,
            closure: self.closure,
            samples: Vec::new(),
            result: None,
        })
    }
    /// Discard an appended suffix and invalidate any cached result. Native
    /// owners use a saved sample count to roll back an entire gauge bank when
    /// a later gauge refuses; ordinary steps need not clone earlier history.
    /// This is an explicit history edit, not a new cyclic loading assumption.
    pub fn discard_after(&mut self, sample_count: u64) -> Result<(), FindleyError> {
        let count = usize::try_from(sample_count).map_err(|_| FindleyError::Budget)?;
        if count > self.samples.len() {
            return Err(FindleyError::Invalid("cannot restore an unobserved prefix"));
        }
        self.samples.truncate(count);
        self.result = None;
        Ok(())
    }
    pub fn push(&mut self, stress_pa: [f64; 6]) -> Result<(), FindleyError> {
        self.push_with_cancel(stress_pa, &mut || false)
    }
    pub fn push_with_cancel(
        &mut self,
        stress_pa: [f64; 6],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FindleyError> {
        self.push_block(core::slice::from_ref(&stress_pa), cancelled)
    }
    pub fn push_block(
        &mut self,
        stress_pa: &[[f64; 6]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FindleyError> {
        if cancelled() {
            return Err(FindleyError::Cancelled);
        }
        if self.result.is_some() {
            return Err(FindleyError::Finished);
        }
        self.admit_count(stress_pa.len())?;
        let previous = self.samples.len();
        for (offset, sample) in stress_pa.iter().enumerate() {
            if cancelled() {
                return Err(FindleyError::Cancelled);
            }
            for (component, value) in sample.iter().enumerate() {
                if !value.is_finite() {
                    return Err(FindleyError::InvalidSample {
                        sample: (previous + offset) as u64,
                        component,
                    });
                }
            }
        }
        self.samples
            .try_reserve(stress_pa.len())
            .map_err(|_| FindleyError::Allocation)?;
        // Roll back only appended values on cancellation; accepted history is
        // never cloned merely to append a bounded block.
        for chunk in stress_pa.chunks(64) {
            if cancelled() {
                self.samples.truncate(previous);
                return Err(FindleyError::Cancelled);
            }
            self.samples.extend_from_slice(chunk);
        }
        if cancelled() {
            self.samples.truncate(previous);
            return Err(FindleyError::Cancelled);
        }
        Ok(())
    }
    fn admit_count(&self, additional: usize) -> Result<(), FindleyError> {
        let next = self
            .samples
            .len()
            .checked_add(additional)
            .ok_or(FindleyError::Budget)?;
        let count = u64::try_from(next).map_err(|_| FindleyError::Budget)?;
        let planes = u64::try_from(self.planes.len()).map_err(|_| FindleyError::Budget)?;
        let work = count.checked_mul(planes).ok_or(FindleyError::Budget)?;
        if count > self.limits.maximum_samples || work > self.limits.maximum_projection_work {
            return Err(FindleyError::Budget);
        }
        Ok(())
    }
    pub fn finish(&mut self) -> Result<FindleyCycleResult, FindleyError> {
        self.finish_with_cancel(&mut || false)
    }
    pub fn finish_with_cancel(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FindleyCycleResult, FindleyError> {
        if cancelled() {
            return Err(FindleyError::Cancelled);
        }
        if let Some(result) = self.result {
            return Ok(result);
        }
        if self.samples.len() < 3 {
            return Err(FindleyError::Invalid(
                "a closed block requires at least three observations including both endpoints",
            ));
        }
        let closure_defect = stress_closure_defect(
            self.samples[0],
            self.samples[self.samples.len() - 1],
            self.closure,
        )?;
        let (origin, scale) = tensor_normalization(&self.samples, cancelled)?;
        let mut path = Vec::new();
        path.try_reserve_exact(self.samples.len())
            .map_err(|_| FindleyError::Allocation)?;
        let mut work = CircleWork {
            used: 0,
            maximum: self.limits.maximum_circle_work,
            cancelled,
        };
        let mut selected: Option<FindleyCycleResult> = None;
        for (plane_index, plane) in self.planes.iter().enumerate() {
            path.clear();
            let mut minimum = f64::INFINITY;
            let mut maximum = f64::NEG_INFINITY;
            for &sample in &self.samples {
                work.check()?;
                let stress = normal_stress(plane.normal_weights, sample)?;
                minimum = minimum.min(stress);
                maximum = maximum.max(stress);
                // Translation in tensor space is a common translation of the
                // shear path, removing static shear before the circle solve.
                let shifted = core::array::from_fn(|i| (sample[i] - origin[i]) / scale);
                let point = [
                    dot6(plane.shear_weights[0], shifted),
                    dot6(plane.shear_weights[1], shifted),
                ];
                if point.iter().any(|v| !v.is_finite()) {
                    return Err(FatigueError::Arithmetic.into());
                }
                path.push(point);
            }
            let normalized_amplitude = enclosing_radius(&mut path, &mut work)?;
            let amplitude = normalized_amplitude * scale;
            let raw = amplitude + self.calibration.normal_stress_sensitivity * maximum;
            let parameter = raw.max(0.0);
            let usage = parameter / self.calibration.limit_factor_pa;
            if !amplitude.is_finite()
                || !raw.is_finite()
                || !usage.is_finite()
                || (normalized_amplitude > 0.0 && amplitude == 0.0)
                || (parameter > 0.0 && usage == 0.0)
            {
                return Err(FatigueError::Arithmetic.into());
            }
            // Select the signed criterion first, then apply the named clipping
            // convention. Exact ties retain the earliest declared plane.
            if selected
                .as_ref()
                .is_none_or(|r| raw > r.unclipped_findley_parameter_pa)
            {
                selected = Some(FindleyCycleResult {
                    samples: self.sample_count(),
                    evaluated_planes: self.planes.len(),
                    selected_plane: plane_index,
                    normal_in_frame: plane.normal,
                    shear_amplitude_pa: amplitude,
                    minimum_normal_stress_pa: minimum,
                    maximum_normal_stress_pa: maximum,
                    unclipped_findley_parameter_pa: raw,
                    findley_parameter_pa: parameter,
                    usage_factor: usage,
                    projection_work: self.projection_work(),
                    circle_work: 0,
                    closure: self.closure,
                    maximum_stress_closure_defect_pa: closure_defect,
                });
            }
        }
        work.check()?;
        let mut result = selected.ok_or(FindleyError::Invalid("missing search planes"))?;
        result.circle_work = work.used;
        self.result = Some(result);
        Ok(result)
    }
}

fn stress_closure_defect(
    first: [f64; 6],
    last: [f64; 6],
    tolerance: StressClosureTolerance,
) -> Result<f64, FindleyError> {
    let mut defect = 0.0_f64;
    for component in 0..6 {
        let (a, b) = (first[component], last[component]);
        if a == b {
            continue;
        }
        let difference = (a - b).abs();
        if !difference.is_finite() {
            return Err(FatigueError::Arithmetic.into());
        }
        let magnitude = a.abs().max(b.abs());
        let scale = magnitude.max(tolerance.absolute_stress_pa);
        let threshold =
            tolerance.absolute_stress_pa / scale + tolerance.relative * (magnitude / scale);
        if difference / scale > threshold {
            return Err(FindleyError::NotClosed { component });
        }
        defect = defect.max(difference);
    }
    Ok(defect)
}
fn midpoint(a: f64, b: f64) -> f64 {
    let span = b - a;
    if span.is_finite() {
        a + span * 0.5
    } else {
        a * 0.5 + b * 0.5
    }
}
fn tensor_normalization(
    samples: &[[f64; 6]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<([f64; 6], f64), FindleyError> {
    let mut lo = samples[0];
    let mut hi = samples[0];
    for sample in samples {
        if cancelled() {
            return Err(FindleyError::Cancelled);
        }
        for i in 0..6 {
            lo[i] = lo[i].min(sample[i]);
            hi[i] = hi[i].max(sample[i]);
        }
    }
    if lo == hi {
        return Err(FindleyError::NotCyclic);
    }
    let origin = core::array::from_fn(|i| midpoint(lo[i], hi[i]));
    let mut scale = 0.0_f64;
    for i in 0..6 {
        scale = scale
            .max((hi[i] - origin[i]).abs())
            .max((lo[i] - origin[i]).abs());
    }
    if !scale.is_finite() || scale <= 0.0 {
        return Err(FatigueError::Arithmetic.into());
    }
    Ok((origin, scale))
}

struct CircleWork<'a, F> {
    used: u64,
    maximum: u64,
    cancelled: &'a mut F,
}
impl<F: FnMut() -> bool> CircleWork<'_, F> {
    fn check(&mut self) -> Result<(), FindleyError> {
        if (self.cancelled)() {
            Err(FindleyError::Cancelled)
        } else {
            Ok(())
        }
    }
    fn step(&mut self) -> Result<(), FindleyError> {
        self.check()?;
        self.used = self.used.checked_add(1).ok_or(FindleyError::Budget)?;
        if self.used > self.maximum {
            Err(FindleyError::Budget)
        } else {
            Ok(())
        }
    }
}
#[derive(Clone, Copy)]
struct Circle {
    center: [f64; 2],
    radius_squared: f64,
}
fn distance_squared(a: [f64; 2], b: [f64; 2]) -> f64 {
    let x = a[0] - b[0];
    let y = a[1] - b[1];
    x * x + y * y
}
fn cross2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
fn subtract2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}
impl Circle {
    fn contains(
        self,
        p: [f64; 2],
        work: &mut CircleWork<'_, impl FnMut() -> bool>,
    ) -> Result<bool, FindleyError> {
        work.step()?;
        Ok(distance_squared(self.center, p) <= self.radius_squared * (1.0 + 32.0 * f64::EPSILON))
    }
    fn diameter(
        a: [f64; 2],
        b: [f64; 2],
        work: &mut CircleWork<'_, impl FnMut() -> bool>,
    ) -> Result<Self, FindleyError> {
        work.step()?;
        let center = [midpoint(a[0], b[0]), midpoint(a[1], b[1])];
        Ok(Self {
            center,
            radius_squared: distance_squared(center, a).max(distance_squared(center, b)),
        })
    }
    fn through_three(
        a: [f64; 2],
        b: [f64; 2],
        c: [f64; 2],
        work: &mut CircleWork<'_, impl FnMut() -> bool>,
    ) -> Result<Self, FindleyError> {
        work.step()?;
        let u = subtract2(b, a);
        let v = subtract2(c, a);
        let determinant = 2.0 * cross2(u, v);
        if determinant == 0.0 {
            return Err(FindleyError::NumericalCircle);
        }
        let uu = u[0] * u[0] + u[1] * u[1];
        let vv = v[0] * v[0] + v[1] * v[1];
        let center = [
            a[0] + (v[1] * uu - u[1] * vv) / determinant,
            a[1] + (u[0] * vv - v[0] * uu) / determinant,
        ];
        let radius_squared = distance_squared(center, a)
            .max(distance_squared(center, b))
            .max(distance_squared(center, c));
        if !radius_squared.is_finite() || center.iter().any(|v| !v.is_finite()) {
            return Err(FindleyError::NumericalCircle);
        }
        Ok(Self {
            center,
            radius_squared,
        })
    }
}

// With a and b constrained to the boundary, the admissible circle centers lie
// on their perpendicular bisector. Each other point imposes a one-sided bound;
// keep the most restrictive circle center on each side, then choose the smaller
// feasible circle. The incremental outer solve supplies the boundary invariant.
fn circle_with_pair(
    points: &[[f64; 2]],
    a: [f64; 2],
    b: [f64; 2],
    work: &mut CircleWork<'_, impl FnMut() -> bool>,
) -> Result<Circle, FindleyError> {
    let diameter = Circle::diameter(a, b, work)?;
    let chord = subtract2(b, a);
    let mut left: Option<(f64, Circle)> = None;
    let mut right: Option<(f64, Circle)> = None;
    for &point in points {
        if diameter.contains(point, work)? {
            continue;
        }
        let side = cross2(chord, subtract2(point, a));
        let circle = Circle::through_three(a, b, point, work)?;
        let location = cross2(chord, subtract2(circle.center, a));
        if side > 0.0 && left.as_ref().is_none_or(|(bound, _)| location > *bound) {
            left = Some((location, circle));
        } else if side < 0.0 && right.as_ref().is_none_or(|(bound, _)| location < *bound) {
            right = Some((location, circle));
        }
    }
    Ok(match (left, right) {
        (None, None) => diameter,
        (Some((_, circle)), None) | (None, Some((_, circle))) => circle,
        (Some((_, left)), Some((_, right))) => {
            if left.radius_squared <= right.radius_squared {
                left
            } else {
                right
            }
        }
    })
}

fn enclosing_radius(
    points: &mut [[f64; 2]],
    work: &mut CircleWork<'_, impl FnMut() -> bool>,
) -> Result<f64, FindleyError> {
    let mut lo = points[0];
    let mut hi = points[0];
    for point in points.iter() {
        work.step()?;
        for i in 0..2 {
            lo[i] = lo[i].min(point[i]);
            hi[i] = hi[i].max(point[i]);
        }
    }
    let origin = [midpoint(lo[0], hi[0]), midpoint(lo[1], hi[1])];
    let scale = (hi[0] - origin[0])
        .abs()
        .max((lo[0] - origin[0]).abs())
        .max((hi[1] - origin[1]).abs())
        .max((lo[1] - origin[1]).abs());
    if scale == 0.0 {
        return Ok(0.0);
    }
    for point in points.iter_mut() {
        work.step()?;
        point[0] = (point[0] - origin[0]) / scale;
        point[1] = (point[1] - origin[1]) / scale;
    }
    // Fixed local permutation improves ordinary incremental workloads. It is
    // deterministic, not a stochastic accuracy claim or a worst-case promise.
    let mut state = 0x6a09_e667_f3bc_c909_u64;
    for i in (1..points.len()).rev() {
        work.step()?;
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let j = (state % (i as u64 + 1)) as usize;
        points.swap(i, j);
    }
    let mut circle = Circle {
        center: points[0],
        radius_squared: 0.0,
    };
    for i in 1..points.len() {
        if circle.contains(points[i], work)? {
            continue;
        }
        work.step()?;
        circle = Circle {
            center: points[i],
            radius_squared: 0.0,
        };
        for j in 0..i {
            if circle.contains(points[j], work)? {
                continue;
            }
            circle = if circle.radius_squared == 0.0 {
                Circle::diameter(points[i], points[j], work)?
            } else {
                circle_with_pair(&points[..j], points[i], points[j], work)?
            };
        }
    }
    let mut enclosing_squared = circle.radius_squared;
    for &point in points.iter() {
        work.step()?;
        let d = distance_squared(circle.center, point);
        if !d.is_finite() || d > circle.radius_squared * (1.0 + 128.0 * f64::EPSILON) {
            return Err(FindleyError::NumericalCircle);
        }
        enclosing_squared = enclosing_squared.max(d);
    }
    // The normalized bounding box is inside a radius-sqrt(2) circle. A larger
    // result is a numerical failure, not an acceptable non-minimal enclosure.
    if enclosing_squared > 2.0 * (1.0 + 128.0 * f64::EPSILON) {
        return Err(FindleyError::NumericalCircle);
    }
    let mut radius = det::sqrt(enclosing_squared);
    if radius * radius < enclosing_squared {
        radius = f64::from_bits(radius.to_bits() + 1);
    }
    let result = radius * scale;
    if !result.is_finite() || (radius > 0.0 && result == 0.0) {
        return Err(FindleyError::NumericalCircle);
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scope() -> FatigueScope {
        FatigueScope {
            source_reference: "synthetic analytic Findley fixture; not measured".into(),
            material: "synthetic elastic coupon".into(),
            process: "fixture".into(),
            surface_and_geometry: "smooth coupon".into(),
            environment: "fixture".into(),
            load_basis: "resolved periodic tensor stress; Findley circle amplitude".into(),
            temperature_k: 293.15,
        }
    }
    fn limits() -> FindleyLimits {
        FindleyLimits {
            maximum_samples: 10_000,
            maximum_planes: 20_000,
            maximum_projection_work: 10_000_000,
            maximum_circle_work: 10_000_000,
        }
    }
    fn cycle(k: f64, f: f64, normals: Vec<[f64; 3]>) -> FindleyCycle {
        FindleyCycle::new(
            Arc::new(FindleyCalibration::new(scope(), k, f).unwrap()),
            normals,
            limits(),
            StressClosureTolerance {
                absolute_stress_pa: 0.0,
                relative: 0.0,
            },
        )
        .unwrap()
    }
    fn evaluate(
        k: f64,
        f: f64,
        normals: Vec<[f64; 3]>,
        samples: &[[f64; 6]],
    ) -> FindleyCycleResult {
        let mut model = cycle(k, f, normals);
        model.push_block(samples, &mut || false).unwrap();
        model.finish().unwrap()
    }
    fn shear(a: f64, b: f64, normal: f64) -> [f64; 6] {
        [normal, 0.0, 0.0, a, 0.0, b]
    }
    fn close(a: f64, b: f64, relative: f64) {
        assert!(
            (a - b).abs() <= relative * a.abs().max(b.abs()).max(f64::MIN_POSITIVE),
            "{a:?} != {b:?}"
        );
    }
    fn circle_history(radius: f64, offset: [f64; 2]) -> Vec<[f64; 6]> {
        let mut samples = Vec::new();
        for i in 0..16 {
            let (s, c) = det::sin_cos(i as f64 * core::f64::consts::TAU / 16.0);
            samples.push(shear(offset[0] + radius * c, offset[1] + radius * s, 0.0));
        }
        samples.push(samples[0]);
        samples
    }

    #[test]
    fn reversing_line_and_separate_normal_peak_use_the_entire_block() {
        let samples = [
            shear(-80.0, 0.0, 10.0),
            shear(0.0, 0.0, 50.0),
            shear(80.0, 0.0, 10.0),
            shear(0.0, 0.0, 10.0),
            shear(-80.0, 0.0, 10.0),
        ];
        let result = evaluate(0.2, 100.0, vec![[1.0, 0.0, 0.0]], &samples);
        close(result.shear_amplitude_pa, 80.0, 2e-14);
        close(result.maximum_normal_stress_pa, 50.0, 2e-14);
        close(result.findley_parameter_pa, 90.0, 2e-14);
        close(result.usage_factor, 0.9, 2e-14);
        assert_eq!(result.samples, 5);
        assert_eq!(result.projection_work, 5);
    }

    #[test]
    fn rotating_shear_with_constant_magnitude_is_not_zero_amplitude() {
        let result = evaluate(
            0.0,
            200.0,
            vec![[1.0, 0.0, 0.0]],
            &circle_history(100.0, [0.0, 0.0]),
        );
        close(result.shear_amplitude_pa, 100.0, 2e-14);
        close(result.usage_factor, 0.5, 2e-14);
    }

    #[test]
    fn triangle_uses_circumcircle_not_half_longest_chord() {
        let height = 50.0 * det::sqrt(3.0);
        let samples = [
            shear(100.0, 0.0, 0.0),
            shear(-50.0, height, 0.0),
            shear(-50.0, -height, 0.0),
            shear(100.0, 0.0, 0.0),
        ];
        let result = evaluate(0.0, 100.0, vec![[1.0, 0.0, 0.0]], &samples);
        close(result.shear_amplitude_pa, 100.0, 3e-14);
        let half_longest_chord = 50.0 * det::sqrt(3.0);
        assert!(result.shear_amplitude_pa > 1.15 * half_longest_chord);
    }

    #[test]
    fn static_shear_translation_and_wide_pressure_scales_preserve_amplitude() {
        assert_eq!(
            normal_stress(
                [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [1e-200, 0.0, 0.0, 1e200, 0.0, 0.0]
            )
            .unwrap(),
            1e-200,
        );
        let original = evaluate(
            0.0,
            200.0,
            vec![[1.0, 0.0, 0.0]],
            &circle_history(100.0, [0.0, 0.0]),
        );
        let shifted = evaluate(
            0.0,
            200.0,
            vec![[1.0, 0.0, 0.0]],
            &circle_history(100.0, [4096.0, -8192.0]),
        );
        close(
            original.shear_amplitude_pa,
            shifted.shear_amplitude_pa,
            3e-13,
        );
        for radius in [1e-200, 1e200] {
            let result = evaluate(
                0.0,
                radius,
                vec![[1.0, 0.0, 0.0]],
                &circle_history(radius, [0.0, 0.0]),
            );
            close(result.shear_amplitude_pa, radius, 3e-14);
            close(result.usage_factor, 1.0, 3e-14);
        }
    }

    #[test]
    fn full_compression_uses_named_zero_clipping_but_static_history_refuses() {
        let samples = [
            shear(-20.0, 0.0, -1000.0),
            shear(20.0, 0.0, -1000.0),
            shear(-20.0, 0.0, -1000.0),
        ];
        let result = evaluate(0.2, 100.0, vec![[1.0, 0.0, 0.0]], &samples);
        close(result.unclipped_findley_parameter_pa, -180.0, 2e-14);
        assert_eq!(result.findley_parameter_pa, 0.0);
        assert_eq!(result.usage_factor, 0.0);
        let mut model = cycle(0.2, 100.0, vec![[1.0, 0.0, 0.0]]);
        model
            .push_block(&[shear(20.0, 30.0, 40.0); 3], &mut || false)
            .unwrap();
        assert_eq!(model.finish(), Err(FindleyError::NotCyclic));
        assert!(model.result().is_none());
    }

    #[test]
    fn uniaxial_and_torsion_match_continuous_analytic_critical_planes() {
        let k = 0.2;
        let root = det::sqrt(1.0 + k * k);
        let axial_plane = [
            det::sqrt((1.0 + k / root) / 2.0),
            det::sqrt((1.0 - k / root) / 2.0),
            0.0,
        ];
        let axial = evaluate(
            k,
            100.0,
            vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], axial_plane],
            &[
                [-100.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [100.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [-100.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            ],
        );
        assert_eq!(axial.selected_plane, 2);
        close(axial.findley_parameter_pa, 50.0 * (root + k), 5e-14);
        let torsion_plane = [
            det::sqrt((1.0 + 1.0 / root) / 2.0),
            det::sqrt((1.0 - 1.0 / root) / 2.0),
            0.0,
        ];
        let torsion = evaluate(
            k,
            100.0,
            vec![[1.0, 0.0, 0.0], torsion_plane],
            &[
                shear(-100.0, 0.0, 0.0),
                shear(100.0, 0.0, 0.0),
                shear(-100.0, 0.0, 0.0),
            ],
        );
        assert_eq!(torsion.selected_plane, 1);
        close(torsion.findley_parameter_pa, 100.0 * root, 5e-14);
    }

    #[test]
    fn shared_hemisphere_refinement_approaches_analytic_plane() {
        let samples = [
            [-100.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [100.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [-100.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ];
        let mut previous = 0.0;
        let analytic = 50.0 * (det::sqrt(1.04) + 0.2);
        for resolution in [2, 4, 8] {
            let planes =
                crate::fatigue::critical_plane::hemisphere_planes(resolution, 1000).unwrap();
            let result = evaluate(0.2, 100.0, planes, &samples);
            assert!(result.findley_parameter_pa + 1e-12 >= previous);
            assert!(result.findley_parameter_pa <= analytic + 1e-12);
            previous = result.findley_parameter_pa;
        }
        assert!(previous > 0.99 * analytic);
    }

    fn rotate_tensor(a: [f64; 6], r: [[f64; 3]; 3]) -> [f64; 6] {
        let matrix = [[a[0], a[3], a[5]], [a[3], a[1], a[4]], [a[5], a[4], a[2]]];
        let mut rotated = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        rotated[i][j] += r[i][k] * matrix[k][l] * r[j][l];
                    }
                }
            }
        }
        [
            rotated[0][0],
            rotated[1][1],
            rotated[2][2],
            rotated[0][1],
            rotated[1][2],
            rotated[2][0],
        ]
    }
    #[test]
    fn common_frame_rotation_preserves_usage_and_selected_plane() {
        let (s, c) = det::sin_cos(0.37);
        let rotation = [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]];
        let planes = vec![
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 2.0, 0.0],
        ];
        let rotated_planes = planes
            .iter()
            .map(|n| core::array::from_fn(|i| (0..3).map(|j| rotation[i][j] * n[j]).sum()))
            .collect();
        let mut samples = Vec::new();
        for i in 0..12 {
            let (a, b) = det::sin_cos(i as f64 * core::f64::consts::TAU / 12.0);
            samples.push([80.0 * b, 20.0 * a, -5.0, 40.0 * a, 15.0 * b, 10.0 * a]);
        }
        samples.push(samples[0]);
        let rotated: Vec<_> = samples
            .iter()
            .map(|&a| rotate_tensor(a, rotation))
            .collect();
        let before = evaluate(0.2, 100.0, planes, &samples);
        let after = evaluate(0.2, 100.0, rotated_planes, &rotated);
        assert_eq!(before.selected_plane, after.selected_plane);
        close(before.shear_amplitude_pa, after.shear_amplitude_pa, 2e-13);
        close(before.usage_factor, after.usage_factor, 2e-13);
    }

    #[test]
    fn closure_checks_unprojected_components_and_preserves_declared_defect() {
        let mut samples = [
            shear(-10.0, 0.0, 0.0),
            shear(10.0, 0.0, 0.0),
            shear(-10.0, 0.0, 0.0),
        ];
        samples[2][2] = 1e-6;
        let mut model = cycle(0.2, 100.0, vec![[1.0, 0.0, 0.0]]);
        model.push_block(&samples, &mut || false).unwrap();
        assert_eq!(
            model.finish(),
            Err(FindleyError::NotClosed { component: 2 })
        );
        assert!(model.result().is_none());
        model.closure.absolute_stress_pa = 1e-6;
        let result = model.finish().unwrap();
        assert_eq!(result.maximum_stress_closure_defect_pa, 1e-6);
        close(result.shear_amplitude_pa, 10.0, 3e-14);
    }

    #[test]
    fn exact_work_caps_refuse_without_publishing_partial_results() {
        let samples = circle_history(100.0, [0.0, 0.0]);
        let baseline = evaluate(0.2, 100.0, vec![[1.0, 0.0, 0.0]], &samples);
        let mut model = cycle(0.2, 100.0, vec![[1.0, 0.0, 0.0]]);
        model.limits.maximum_samples = samples.len() as u64;
        model.limits.maximum_projection_work = samples.len() as u64;
        model.limits.maximum_circle_work = baseline.circle_work;
        model.push_block(&samples, &mut || false).unwrap();
        assert_eq!(model.push(samples[0]), Err(FindleyError::Budget));
        assert_eq!(model.finish().unwrap(), baseline);
        let mut insufficient = model.fresh().unwrap();
        insufficient.limits.maximum_circle_work -= 1;
        insufficient.push_block(&samples, &mut || false).unwrap();
        assert_eq!(insufficient.finish(), Err(FindleyError::Budget));
        assert_eq!(insufficient.sample_count(), samples.len() as u64);
        assert!(insufficient.result().is_none());
        let mut projections = model.fresh().unwrap();
        projections.limits.maximum_projection_work -= 1;
        assert_eq!(
            projections.push_block(&samples, &mut || false),
            Err(FindleyError::Budget)
        );
        assert_eq!(projections.sample_count(), 0);
    }

    #[test]
    fn cancellation_invalid_samples_and_native_prefix_rollback_are_atomic() {
        let mut model = cycle(0.2, 100.0, vec![[1.0, 0.0, 0.0]]);
        let start = shear(-10.0, 0.0, 0.0);
        model.push(start).unwrap();
        let mut calls = 0;
        assert_eq!(
            model.push_block(&[shear(10.0, 0.0, 0.0); 10], &mut || {
                calls += 1;
                calls == 13
            }),
            Err(FindleyError::Cancelled)
        );
        assert_eq!(model.sample_count(), 1);
        assert_eq!(model.samples[0], start);
        let invalid = [shear(10.0, 0.0, 0.0), shear(f64::NAN, 0.0, 0.0)];
        assert_eq!(
            model.push_block(&invalid, &mut || false),
            Err(FindleyError::InvalidSample {
                sample: 2,
                component: 3
            })
        );
        assert_eq!(model.sample_count(), 1);
        model
            .push_block(&[shear(10.0, 0.0, 0.0), start], &mut || false)
            .unwrap();
        let mut calls = 0;
        assert_eq!(
            model.finish_with_cancel(&mut || {
                calls += 1;
                calls == 12
            }),
            Err(FindleyError::Cancelled)
        );
        assert!(model.result().is_none());
        assert_eq!(model.sample_count(), 3);
        let original = model.finish().unwrap();
        assert_eq!(model.push(start), Err(FindleyError::Finished));
        assert!(model.discard_after(4).is_err());
        assert_eq!(model.result(), Some(&original));
        model.discard_after(2).unwrap();
        assert!(model.result().is_none());
        model.push(start).unwrap();
        assert_eq!(model.finish().unwrap(), original);
        let fresh = model.fresh().unwrap();
        assert_eq!(fresh.sample_count(), 0);
        assert!(fresh.result().is_none());
        assert!(Arc::ptr_eq(&fresh.calibration, &model.calibration));
    }

    // Independent tiny-set oracle: enumerate all two/three-point support
    // circles, retain only enclosures, then minimize radius. Unlike the kernel,
    // this does no incremental boundary updates or fixed permutation.
    fn support_oracle(points: &[[f64; 2]]) -> f64 {
        let mut candidates = Vec::new();
        for i in 0..points.len() {
            for j in 0..i {
                let center = [
                    (points[i][0] + points[j][0]) / 2.0,
                    (points[i][1] + points[j][1]) / 2.0,
                ];
                candidates.push((center, distance_squared(center, points[i])));
                for k in 0..j {
                    let [ax, ay] = points[i];
                    let [bx, by] = points[j];
                    let [cx, cy] = points[k];
                    let denominator = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
                    if denominator == 0.0 {
                        continue;
                    }
                    let aa = ax * ax + ay * ay;
                    let bb = bx * bx + by * by;
                    let cc = cx * cx + cy * cy;
                    let center = [
                        (aa * (by - cy) + bb * (cy - ay) + cc * (ay - by)) / denominator,
                        (aa * (cx - bx) + bb * (ax - cx) + cc * (bx - ax)) / denominator,
                    ];
                    candidates.push((center, distance_squared(center, points[i])));
                }
            }
        }
        let best = candidates
            .into_iter()
            .filter(|(center, r2)| {
                points
                    .iter()
                    .all(|p| distance_squared(*center, *p) <= r2 * (1.0 + 1e-12))
            })
            .map(|(_, r2)| r2)
            .fold(f64::INFINITY, f64::min);
        det::sqrt(best)
    }
    #[test]
    fn enclosing_circle_matches_independent_support_enumeration() {
        let mut cases = vec![
            vec![[-2.0, 0.0], [0.0, 0.0], [3.0, 0.0], [1.0, 0.0]],
            vec![[-2.0, 0.0], [0.0, 1e-12], [3.0, 0.0], [1.0, -1e-12]],
            vec![[0.0, 0.0], [4.0, 0.0], [1.0, 0.2], [0.0, 0.0]],
        ];
        for case in 0..24 {
            cases.push(
                (0..7)
                    .map(|i| {
                        let (s, c) = det::sin_cos((i * 7 + case * 13) as f64 * 0.73);
                        [c * (1.0 + 0.1 * i as f64), s * (1.0 + 0.03 * case as f64)]
                    })
                    .collect(),
            );
        }
        for points in cases {
            let expected = support_oracle(&points);
            let mut normalized = points.clone();
            let mut cancelled = || false;
            let mut work = CircleWork {
                used: 0,
                maximum: 100_000,
                cancelled: &mut cancelled,
            };
            let actual = enclosing_radius(&mut normalized, &mut work).unwrap();
            close(actual, expected, 2e-12);
        }
    }
}

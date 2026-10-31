//! Private fixed-plane projections and bounded two-dimensional shear paths.
//! Shared by stress-based Findley and strain-based Fatemi-Socie; model-specific
//! calibration, plane selection and normal-stress treatment stay with each law.

use crate::fatigue::{FatigueError, NormalStressDirection};
use fs_math::det;

#[derive(Debug, Clone, PartialEq)]
pub(super) enum ShearPathError {
    Budget,
    Cancelled,
    NotCyclic,
    NumericalCircle,
    Fatigue(FatigueError),
}
impl From<FatigueError> for ShearPathError {
    fn from(error: FatigueError) -> Self {
        Self::Fatigue(error)
    }
}

#[derive(Debug, Clone, Copy)]
pub(super) struct SearchPlane {
    pub(super) normal: [f64; 3],
    pub(super) normal_weights: [f64; 6],
    pub(super) shear_weights: [[f64; 6]; 2],
}
impl SearchPlane {
    pub(super) fn new(direction: [f64; 3]) -> Result<Self, ShearPathError> {
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
pub(super) fn dot6(a: [f64; 6], b: [f64; 6]) -> f64 {
    let mut result = 0.0;
    for i in 0..6 {
        result += a[i] * b[i];
    }
    result
}
pub(super) fn normal_stress(weights: [f64; 6], stress: [f64; 6]) -> Result<f64, ShearPathError> {
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

pub(super) fn midpoint(a: f64, b: f64) -> f64 {
    let span = b - a;
    if span.is_finite() {
        a + span * 0.5
    } else {
        a * 0.5 + b * 0.5
    }
}
pub(super) fn tensor_normalization(
    samples: impl IntoIterator<Item = [f64; 6]>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<([f64; 6], f64), ShearPathError> {
    let mut samples = samples.into_iter();
    let first = samples.next().ok_or(FatigueError::Arithmetic)?;
    let mut lo = first;
    let mut hi = first;
    for sample in core::iter::once(first).chain(samples) {
        if cancelled() {
            return Err(ShearPathError::Cancelled);
        }
        for i in 0..6 {
            lo[i] = lo[i].min(sample[i]);
            hi[i] = hi[i].max(sample[i]);
        }
    }
    if lo == hi {
        return Err(ShearPathError::NotCyclic);
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

pub(super) struct CircleWork<'a, F> {
    pub(super) used: u64,
    pub(super) maximum: u64,
    pub(super) cancelled: &'a mut F,
}
impl<F: FnMut() -> bool> CircleWork<'_, F> {
    pub(super) fn check(&mut self) -> Result<(), ShearPathError> {
        if (self.cancelled)() {
            Err(ShearPathError::Cancelled)
        } else {
            Ok(())
        }
    }
    fn step(&mut self) -> Result<(), ShearPathError> {
        self.check()?;
        self.used = self.used.checked_add(1).ok_or(ShearPathError::Budget)?;
        if self.used > self.maximum {
            Err(ShearPathError::Budget)
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
pub(super) fn distance_squared(a: [f64; 2], b: [f64; 2]) -> f64 {
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
    ) -> Result<bool, ShearPathError> {
        work.step()?;
        Ok(distance_squared(self.center, p) <= self.radius_squared * (1.0 + 32.0 * f64::EPSILON))
    }
    fn diameter(
        a: [f64; 2],
        b: [f64; 2],
        work: &mut CircleWork<'_, impl FnMut() -> bool>,
    ) -> Result<Self, ShearPathError> {
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
    ) -> Result<Self, ShearPathError> {
        work.step()?;
        let u = subtract2(b, a);
        let v = subtract2(c, a);
        let determinant = 2.0 * cross2(u, v);
        if determinant == 0.0 {
            return Err(ShearPathError::NumericalCircle);
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
            return Err(ShearPathError::NumericalCircle);
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
) -> Result<Circle, ShearPathError> {
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

pub(super) fn enclosing_radius(
    points: &mut [[f64; 2]],
    work: &mut CircleWork<'_, impl FnMut() -> bool>,
) -> Result<f64, ShearPathError> {
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
            return Err(ShearPathError::NumericalCircle);
        }
        enclosing_squared = enclosing_squared.max(d);
    }
    // The normalized bounding box is inside a radius-sqrt(2) circle. A larger
    // result is a numerical failure, not an acceptable non-minimal enclosure.
    if enclosing_squared > 2.0 * (1.0 + 128.0 * f64::EPSILON) {
        return Err(ShearPathError::NumericalCircle);
    }
    let mut radius = det::sqrt(enclosing_squared);
    if radius * radius < enclosing_squared {
        radius = f64::from_bits(radius.to_bits() + 1);
    }
    let result = radius * scale;
    if !result.is_finite() || (radius > 0.0 && result == 0.0) {
        return Err(ShearPathError::NumericalCircle);
    }
    Ok(result)
}

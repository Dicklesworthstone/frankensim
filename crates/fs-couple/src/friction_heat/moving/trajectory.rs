//! Exact swept residence for accelerated translating planar contact.
//!
//! The supplied placements describe an existing mechanical trajectory. They
//! do not integrate forces or define a spherical/rotating material map. Edge
//! crossings split the quadratic overlap into polynomial intervals; Simpson
//! integration is exact for their area and first time moment, up to floating
//! arithmetic. Endpoint thermal coefficients still require time refinement.

use super::{MovingContactPath, MovingContactPosition};

/// Actual placements at the beginning, midpoint and end of one interval.
/// Each offset is interpolated by its unique quadratic in normalized time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadraticContactPath {
    pub start: MovingContactPosition,
    pub midpoint: MovingContactPosition,
    pub end: MovingContactPosition,
}

/// Invalid or unrepresentable supplied motion/work; no owner is modified.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TrajectoryError {
    Invalid(&'static str),
    Arithmetic,
}
impl core::fmt::Display for TrajectoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "moving friction trajectory: {self:?}")
    }
}
impl std::error::Error for TrajectoryError {}

fn finite(x: f64) -> Result<f64, TrajectoryError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(TrajectoryError::Arithmetic)
    }
}

/// Polynomial in normalized time, c0 + c1*t + c2*t².
#[derive(Clone, Copy, Debug)]
pub(super) struct Curve([f64; 3]);
impl Curve {
    fn linear(start: f64, end: f64) -> Result<Self, TrajectoryError> {
        finite(start)?;
        finite(end)?;
        Ok(Self([start, finite(end - start)?, 0.0]))
    }
    fn through(start: f64, midpoint: f64, end: f64) -> Result<Self, TrajectoryError> {
        finite(start)?;
        finite(midpoint)?;
        finite(end)?;
        let middle = finite(midpoint - start)?;
        let change = finite(end - start)?;
        let quadratic = finite(2.0 * finite(change - finite(2.0 * middle)?)?)?;
        Ok(Self([start, finite(change - quadratic)?, quadratic]))
    }
    pub(super) fn shift(self, amount: f64) -> Result<Self, TrajectoryError> {
        Ok(Self([finite(self.0[0] + amount)?, self.0[1], self.0[2]]))
    }
    fn subtract(self, other: Self) -> Result<Self, TrajectoryError> {
        Ok(Self([
            finite(self.0[0] - other.0[0])?,
            finite(self.0[1] - other.0[1])?,
            finite(self.0[2] - other.0[2])?,
        ]))
    }
    fn at(self, time: f64) -> Result<f64, TrajectoryError> {
        finite(self.0[2].mul_add(time, self.0[1]).mul_add(time, self.0[0]))
    }
    pub(super) fn bounds(self) -> Result<[f64; 2], TrajectoryError> {
        let a = self.at(0.0)?;
        let b = self.at(1.0)?;
        let mut bounds = [a.min(b), a.max(b)];
        if self.0[2] != 0.0 {
            // Divide before multiplying to avoid overflow in 2*c2.
            let vertex = -0.5 * (self.0[1] / self.0[2]);
            if vertex > 0.0 && vertex < 1.0 {
                let value = self.at(vertex)?;
                bounds = [bounds[0].min(value), bounds[1].max(value)];
            }
        }
        Ok(bounds)
    }
    fn crossings(self, times: &mut Vec<f64>) {
        let scale = self.0.into_iter().map(f64::abs).fold(0.0_f64, f64::max);
        if scale == 0.0 {
            return;
        }
        let [c, b, a] = self.0.map(|x| x / scale);
        let mut push = |t: f64| {
            if t > 0.0 && t < 1.0 {
                times.push(t);
            }
        };
        if a == 0.0 {
            if b != 0.0 {
                push(-c / b);
            }
            return;
        }
        let discriminant = b.mul_add(b, -4.0 * a * c);
        if discriminant < 0.0 {
            return;
        }
        let root = discriminant.sqrt();
        // This stable pair avoids subtracting nearly equal positive roots.
        let q = -0.5 * (b + root.copysign(b));
        if q == 0.0 {
            push(-0.5 * (b / a));
        } else {
            push(q / a);
            push(c / q);
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) struct ContactCurves {
    pub bodies: [Curve; 2],
    pub footprint: Curve,
    pub end: MovingContactPosition,
}
impl ContactCurves {
    pub(super) fn linear(path: MovingContactPath) -> Result<Self, TrajectoryError> {
        Ok(Self {
            bodies: [
                Curve::linear(path.start.body_offsets_m[0], path.end.body_offsets_m[0])?,
                Curve::linear(path.start.body_offsets_m[1], path.end.body_offsets_m[1])?,
            ],
            footprint: Curve::linear(path.start.footprint_offset_m, path.end.footprint_offset_m)?,
            end: path.end,
        })
    }
    pub(super) fn quadratic(path: QuadraticContactPath) -> Result<Self, TrajectoryError> {
        Ok(Self {
            bodies: [
                Curve::through(
                    path.start.body_offsets_m[0],
                    path.midpoint.body_offsets_m[0],
                    path.end.body_offsets_m[0],
                )?,
                Curve::through(
                    path.start.body_offsets_m[1],
                    path.midpoint.body_offsets_m[1],
                    path.end.body_offsets_m[1],
                )?,
            ],
            footprint: Curve::through(
                path.start.footprint_offset_m,
                path.midpoint.footprint_offset_m,
                path.end.footprint_offset_m,
            )?,
            end: path.end,
        })
    }
    /// Full-interval coverage, including excursions between the three samples.
    pub(super) fn contained(
        self,
        side: usize,
        footprint: [f64; 2],
        surface: [f64; 2],
    ) -> Result<bool, TrajectoryError> {
        let relative = self.footprint.subtract(self.bodies[side])?.bounds()?;
        Ok(finite(footprint[0] + relative[0])? >= surface[0]
            && finite(footprint[1] + relative[1])? <= surface[1])
    }
}
pub(super) fn edges(span: [f64; 2], offset: Curve) -> Result<[Curve; 2], TrajectoryError> {
    Ok([offset.shift(span[0])?, offset.shift(span[1])?])
}
pub(super) fn bounds(edges: [Curve; 2]) -> Result<[f64; 2], TrajectoryError> {
    Ok([edges[0].bounds()?[0], edges[1].bounds()?[1]])
}

/// Integral of overlap and normalized-time times overlap over [0,1].
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Moments {
    pub area: f64,
    pub time_area: f64,
}
fn overlap(edges: &[Curve; 6], time: f64) -> Result<f64, TrajectoryError> {
    let left = edges[0]
        .at(time)?
        .max(edges[2].at(time)?)
        .max(edges[4].at(time)?);
    let right = edges[1]
        .at(time)?
        .min(edges[3].at(time)?)
        .min(edges[5].at(time)?);
    Ok(finite(right - left)?.max(0.0))
}
pub(super) fn swept_moments(edges: [Curve; 6]) -> Result<Moments, TrajectoryError> {
    // Six edges have at most thirty interior pairwise quadratic crossings.
    // Constant bounded work between the original owner's per-pair checkpoints.
    let mut times = Vec::with_capacity(32);
    times.extend([0.0, 1.0]);
    for i in 0..6 {
        for j in i + 1..6 {
            edges[i].subtract(edges[j])?.crossings(&mut times);
        }
    }
    times.sort_by(f64::total_cmp);
    times.dedup();
    let mut result = Moments {
        area: 0.0,
        time_area: 0.0,
    };
    for interval in times.windows(2) {
        let [a, b] = [interval[0], interval[1]];
        let middle = f64::midpoint(a, b);
        let [fa, fm, fb] = [
            overlap(&edges, a)?,
            overlap(&edges, middle)?,
            overlap(&edges, b)?,
        ];
        // Positive weights avoid cancellation for narrow positive residence.
        let factor = (b - a) / 6.0;
        result.area = finite(result.area + factor * (fa + 4.0 * fm + fb))?;
        result.time_area =
            finite(result.time_area + factor * (a * fa + 4.0 * middle * fm + b * fb))?;
    }
    Ok(result)
}

/// Integrate an affine, nonnegative physical power density over actual residence.
pub(super) fn affine_work(
    area: f64,
    duration: f64,
    power: [f64; 2],
    moments: impl ExactSizeIterator<Item = Moments>,
) -> Result<Vec<f64>, TrajectoryError> {
    if !area.is_finite()
        || area <= 0.0
        || !duration.is_finite()
        || duration <= 0.0
        || power.iter().any(|p| !p.is_finite() || *p < 0.0)
        || moments.len() == 0
    {
        return Err(TrajectoryError::Invalid(
            "positive area/duration and nonnegative finite endpoint powers required",
        ));
    }
    let total = finite(duration * finite(0.5 * power[0] + 0.5 * power[1])?)?;
    if total == 0.0 && power != [0.0; 2] {
        return Err(TrajectoryError::Arithmetic);
    }
    let mut result = Vec::with_capacity(moments.len());
    let mut sum = 0.0;
    let mut area_sum = 0.0;
    let mut moment_sum = 0.0;
    let mut largest = 0usize;
    for moment in moments {
        if !moment.area.is_finite()
            || !moment.time_area.is_finite()
            || moment.area <= 0.0
            || moment.time_area < 0.0
            || moment.time_area > moment.area
        {
            return Err(TrajectoryError::Invalid("invalid swept area/time moment"));
        }
        area_sum = finite(area_sum + moment.area)?;
        moment_sum = finite(moment_sum + moment.time_area)?;
        let weight0 = finite((moment.area - moment.time_area) / area)?;
        let weight1 = finite(moment.time_area / area)?;
        let work = finite(duration * finite(power[0] * weight0 + power[1] * weight1)?)?;
        if result.is_empty() || work > result[largest] {
            largest = result.len();
        }
        result.push(work);
        sum = finite(sum + work)?;
    }
    let tolerance = 2048.0 * f64::EPSILON * (result.len() as f64).sqrt();
    if (area_sum - area).abs() > tolerance * area
        || (moment_sum - 0.5 * area).abs() > tolerance * area
        || (sum - total).abs() > tolerance * total
    {
        return Err(TrajectoryError::Invalid(
            "swept physical power does not conserve footprint work",
        ));
    }
    // Correct only accumulated roundoff, at the largest receipt rather than a
    // nearly empty last patch. No inferred physical loss or body partition.
    result[largest] = finite(result[largest] + finite(total - sum)?)?;
    if result[largest] < 0.0 {
        return Err(TrajectoryError::Arithmetic);
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn close(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 2e-13 * a.abs().max(b.abs()).max(1.0),
            "{a} != {b}"
        );
    }
    fn position(x: f64) -> MovingContactPosition {
        MovingContactPosition {
            body_offsets_m: [0.0; 2],
            footprint_offset_m: x,
        }
    }
    fn sweep(span: [f64; 2], footprint: Curve) -> Moments {
        let fixed = Curve::linear(0.0, 0.0).unwrap();
        let a = edges(span, fixed).unwrap();
        let f = edges([0.0, 1.0], footprint).unwrap();
        swept_moments([a[0], a[1], a[0], a[1], f[0], f[1]]).unwrap()
    }
    #[test]
    fn accelerating_native_segment_power_reaches_the_correct_material() {
        // Constant traction=1, slip v=2t, position=t², duration=1.
        let foot = Curve::through(0.0, 0.25, 1.0).unwrap();
        let moments = [sweep([0.0, 1.0], foot), sweep([1.0, 2.0], foot)];
        close(moments[0].area, 2.0 / 3.0);
        close(moments[1].area, 1.0 / 3.0);
        close(moments[0].time_area, 0.25);
        close(moments[1].time_area, 0.25);
        let work = affine_work(1.0, 1.0, [0.0, 2.0], moments.into_iter()).unwrap();
        close(work[0], 0.5);
        close(work[1], 0.5);
        close(work.iter().sum(), 1.0);
        let uniform = affine_work(1.0, 1.0, [1.0; 2], moments.into_iter()).unwrap();
        close(uniform[0], 2.0 / 3.0);
        close(uniform[1], 1.0 / 3.0);
    }
    #[test]
    fn reversal_and_unsampled_extrema_preserve_full_path_geometry() {
        let path = QuadraticContactPath {
            start: position(0.0),
            midpoint: position(1.0),
            end: position(0.0),
        };
        let curve = ContactCurves::quadratic(path).unwrap();
        let moments = [
            sweep([0.0, 1.0], curve.footprint),
            sweep([1.0, 2.0], curve.footprint),
        ];
        close(moments[0].area, 1.0 / 3.0);
        close(moments[1].area, 2.0 / 3.0);
        close(moments[0].time_area, 1.0 / 6.0);
        close(moments[1].time_area, 1.0 / 3.0);
        assert!(curve.contained(0, [0.0, 1.0], [0.0, 2.0]).unwrap());
        assert!(!curve.contained(0, [0.0, 1.0], [0.0, 1.5]).unwrap());
        // All three supplied placements fit, but the vertex t=.625 does not.
        let missed = ContactCurves::quadratic(QuadraticContactPath {
            start: position(0.0),
            midpoint: position(1.5),
            end: position(1.0),
        })
        .unwrap();
        assert!(!missed.contained(0, [0.0, 1.0], [0.0, 2.5]).unwrap());
    }
    #[test]
    fn linear_paths_common_translation_and_endpoint_power_recover_original_limits() {
        let linear = ContactCurves::linear(MovingContactPath {
            start: position(0.0),
            end: position(1.0),
        })
        .unwrap();
        let quadratic = ContactCurves::quadratic(QuadraticContactPath {
            start: position(0.0),
            midpoint: position(0.5),
            end: position(1.0),
        })
        .unwrap();
        let a = [
            sweep([0.0, 1.0], linear.footprint),
            sweep([1.0, 2.0], linear.footprint),
        ];
        let b = [
            sweep([0.0, 1.0], quadratic.footprint),
            sweep([1.0, 2.0], quadratic.footprint),
        ];
        assert_eq!(a, b);
        close(a[0].area, 0.5);
        close(a[0].time_area, 1.0 / 6.0);
        let rising = affine_work(1.0, 2.0, [0.0, 2.0], a.into_iter()).unwrap();
        close(rising[0], 2.0 / 3.0);
        close(rising[1], 4.0 / 3.0);
        let translated = QuadraticContactPath {
            start: MovingContactPosition {
                body_offsets_m: [10.0; 2],
                footprint_offset_m: 10.0,
            },
            midpoint: MovingContactPosition {
                body_offsets_m: [12.0; 2],
                footprint_offset_m: 12.5,
            },
            end: MovingContactPosition {
                body_offsets_m: [13.0; 2],
                footprint_offset_m: 14.0,
            },
        };
        let moved = ContactCurves::quadratic(translated).unwrap();
        for (i, expected) in a.into_iter().enumerate() {
            let faces = edges([i as f64, i as f64 + 1.0], moved.bodies[0]).unwrap();
            let foot = edges([0.0, 1.0], moved.footprint).unwrap();
            let actual =
                swept_moments([faces[0], faces[1], faces[0], faces[1], foot[0], foot[1]]).unwrap();
            close(actual.area, expected.area);
            close(actual.time_area, expected.time_area);
        }
        assert!(affine_work(1.0, 1.0, [-1.0, 2.0], a.into_iter()).is_err());
        assert!(affine_work(1.0, 1.0, [f64::NAN, 0.0], a.into_iter()).is_err());
        assert!(
            ContactCurves::quadratic(QuadraticContactPath {
                start: position(0.0),
                midpoint: position(f64::INFINITY),
                end: position(0.0),
            })
            .is_err()
        );
    }
}

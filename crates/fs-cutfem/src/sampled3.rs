//! Trilinearly reconstructed sampled geometry for every [`CutSdf3`] consumer.
//!
//! The represented domain is the negative set of the continuous interpolant on
//! the retained floating-point grid knots, NOT an unknown surface that happened
//! to produce the samples. Enclosures cover that interpolant and its one-sided
//! partial derivatives, including both sides of grid seams. They do not certify
//! sampling error, Euclidean distance, topology, or quadrature accuracy.
//!
//! Queries never extrapolate: an uncovered box returns an unbounded enclosure
//! and an uncovered point returns NaN, causing the existing solver to refuse.
//! Large queries use precomputed global bounds rather than unbounded scans.

use crate::{CutSdf3, HeightAxis};
use fs_ivl::Interval;
use std::ops::ControlFlow;

/// Admission limit on retained scalar samples, before any internal allocation.
pub const MAX_SAMPLED3_VALUES: usize = 1_048_576;
/// Maximum grid cells inspected by one enclosure callback.
const LOCAL_CELLS: usize = 512;

/// Invalid sampled input or interrupted construction; no partial field escapes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SampledSdf3Error {
    /// The input cannot define the documented finite, covered interpolant.
    Invalid(&'static str),
    /// The caller stopped construction at a bounded work boundary.
    Cancelled,
}
impl std::fmt::Display for SampledSdf3Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "sampled 3-D field refused: {reason}"),
            Self::Cancelled => f.write_str("sampled 3-D field construction cancelled"),
        }
    }
}
impl std::error::Error for SampledSdf3Error {}

/// Immutable, negative-inside trilinear field on a rectilinear sampled grid.
///
/// Input is a uniform grid. Its knots are retained as `origin + i * spacing`;
/// interpolation uses the actual adjacent knot difference, so rounding at a
/// translated origin neither creates gaps nor changes the support silently.
#[derive(Debug, Clone)]
pub struct SampledSdf3 {
    axes: [Vec<f64>; 3],
    shape: [usize; 3],
    samples: Vec<f64>,
    global: Interval,
    slopes: [Interval; 3],
}

fn point(value: f64) -> Interval { Interval::new(value, value) }
fn whole() -> Interval { Interval::new(f64::NEG_INFINITY, f64::INFINITY) }
fn hull(a: Interval, b: Interval) -> Interval {
    Interval::new(a.lo().min(b.lo()), a.hi().max(b.hi()))
}
fn lerp(a: Interval, b: Interval, t: Interval) -> Interval {
    (point(1.0) - t) * a + t * b
}
fn poll(callback: &mut dyn FnMut() -> ControlFlow<()>) -> Result<(), SampledSdf3Error> {
    match callback() {
        ControlFlow::Continue(()) => Ok(()),
        ControlFlow::Break(()) => Err(SampledSdf3Error::Cancelled),
    }
}

impl SampledSdf3 {
    /// Admit samples in x-fastest order: `x + nx * (y + ny * z)`.
    ///
    /// Each axis needs at least two knots, positive spacing in [1e-100, 1e100],
    /// and finite coordinates of magnitude at most 1e100. Samples obey the same
    /// magnitude bound; missing/nonfinite values are never filled or clamped.
    /// Counts, products and distinct rounded knots are checked before use.
    /// The cancellation callback is polled before allocation, every 256 knots
    /// or samples, and before publication. An edge-slope sweep shares those
    /// sample checkpoints; no separate uncancellable full-grid pass is needed.
    pub fn with_control(
        origin: [f64; 3], spacing: [f64; 3], shape: [usize; 3], samples: Vec<f64>,
        callback: &mut dyn FnMut() -> ControlFlow<()>,
    ) -> Result<Self, SampledSdf3Error> {
        poll(callback)?;
        let invalid = SampledSdf3Error::Invalid;
        let count = shape.iter().try_fold(1usize, |n, &d| {
            if d < 2 { None } else { n.checked_mul(d) }
        }).ok_or(invalid("each axis needs at least two knots and a representable sample count"))?;
        if count > MAX_SAMPLED3_VALUES || samples.len() != count {
            return Err(invalid("sample count must match shape and not exceed 1048576"));
        }
        let mut axes: [Vec<f64>; 3] = std::array::from_fn(|_| Vec::new());
        for axis in 0..3 {
            if !origin[axis].is_finite() || origin[axis].abs() > 1e100
                || !spacing[axis].is_finite() || !(1e-100..=1e100).contains(&spacing[axis]) {
                return Err(invalid("origin/spacing is outside the admitted finite range"));
            }
            axes[axis].reserve(shape[axis]);
            for i in 0..shape[axis] {
                if i % 256 == 0 { poll(callback)?; }
                let knot = origin[axis] + (i as f64) * spacing[axis];
                if !knot.is_finite() || knot.abs() > 1e100
                    || axes[axis].last().is_some_and(|&previous| knot <= previous) {
                    return Err(invalid("grid knots must be finite, bounded and strictly increasing"));
                }
                axes[axis].push(knot);
            }
        }
        let mut global = point(0.0);
        let mut slopes = [point(0.0); 3];
        let mut have_slope = [false; 3];
        let strides = [1, shape[0], shape[0] * shape[1]];
        for (at, &sample) in samples.iter().enumerate() {
            if at % 256 == 0 { poll(callback)?; }
            if !sample.is_finite() || sample.abs() > 1e100 {
                return Err(invalid("samples must be finite with magnitude at most 1e100"));
            }
            global = if at == 0 { point(sample) } else { hull(global, point(sample)) };
            let index = [at % shape[0], (at / shape[0]) % shape[1], at / strides[2]];
            for axis in 0..3 {
                if index[axis] == 0 { continue; }
                // Only already-validated preceding samples are read here.
                let width = point(axes[axis][index[axis]]) - point(axes[axis][index[axis] - 1]);
                if width.lo() <= 0.0 { return Err(invalid("adjacent grid spacing is not interval-positive")); }
                let slope = (point(sample) - point(samples[at - strides[axis]])) / width;
                if !slope.lo().is_finite() || !slope.hi().is_finite() {
                    return Err(invalid("an edge derivative cannot be represented finitely"));
                }
                slopes[axis] = if have_slope[axis] { hull(slopes[axis], slope) } else { slope };
                have_slope[axis] = true;
            }
        }
        poll(callback)?;
        Ok(Self { axes, shape, samples, global, slopes })
    }

    /// Exact covered bounds; outside queries do not reuse boundary samples.
    #[must_use]
    pub fn bounds(&self) -> ([f64; 3], [f64; 3]) {
        (std::array::from_fn(|a| self.axes[a][0]),
         std::array::from_fn(|a| self.axes[a][self.shape[a] - 1]))
    }

    fn covered(&self, lo: [f64; 3], hi: [f64; 3]) -> bool {
        (0..3).all(|a| lo[a].is_finite() && hi[a].is_finite() && lo[a] <= hi[a]
            && lo[a] >= self.axes[a][0] && hi[a] <= self.axes[a][self.shape[a] - 1])
    }

    fn offset(&self, index: [usize; 3]) -> usize {
        index[0] + self.shape[0] * (index[1] + self.shape[1] * index[2])
    }

    fn corners(&self, cell: [usize; 3]) -> [f64; 8] {
        std::array::from_fn(|corner| {
            self.samples[self.offset(std::array::from_fn(|a| cell[a] + ((corner >> a) & 1)))]
        })
    }

    fn cell_range(&self, cell: [usize; 3], p: [f64; 3], derivative: Option<usize>) -> Interval {
        let mut values = self.corners(cell).map(point);
        for a in 0..3 {
            let left = point(self.axes[a][cell[a]]);
            let width = point(self.axes[a][cell[a] + 1]) - left;
            let t = (point(p[a]) - left) / width;
            // Membership of p in this closed cell independently proves [0,1].
            let t = Interval::new(t.lo().max(0.0), t.hi().min(1.0));
            let stride = 1 << a;
            for at in (0..8).step_by(2 * stride) {
                values[at] = if derivative == Some(a) {
                    (values[at + stride] - values[at]) / width
                } else { lerp(values[at], values[at + stride], t) };
            }
        }
        values[0]
    }

    fn range(&self, lo: [f64; 3], hi: [f64; 3], derivative: Option<usize>) -> Interval {
        if !self.covered(lo, hi) { return whole(); }
        // Include BOTH adjacent cells for exact-knot endpoints. In particular,
        // a degenerate point box on a seam must retain both one-sided slopes.
        let first: [usize; 3] = std::array::from_fn(|a|
            self.axes[a].partition_point(|&x| x < lo[a]).saturating_sub(1).min(self.shape[a] - 2));
        let last: [usize; 3] = std::array::from_fn(|a|
            self.axes[a].partition_point(|&x| x <= hi[a]).saturating_sub(1).min(self.shape[a] - 2));
        let global = derivative.map_or(self.global, |a| self.slopes[a]);
        let cells: usize = (0..3).map(|a| last[a] - first[a] + 1).product();
        if cells > LOCAL_CELLS { return global; }
        let mut result = None;
        for z in first[2]..=last[2] {
            for y in first[1]..=last[1] {
                for x in first[0]..=last[0] {
                    let cell = [x, y, z];
                    let left: [f64; 3] = std::array::from_fn(|a| lo[a].max(self.axes[a][cell[a]]));
                    let right: [f64; 3] = std::array::from_fn(|a| hi[a].min(self.axes[a][cell[a] + 1]));
                    // A multi-affine function (or its bilinear partial) attains
                    // its extrema at vertices of each clipped Cartesian box.
                    for corner in 0..8 {
                        if (0..3).any(|a| corner & (1 << a) != 0
                            && (left[a] == right[a] || derivative == Some(a))) { continue; }
                        let p = std::array::from_fn(|a| if corner & (1 << a) == 0 { left[a] } else { right[a] });
                        let value = self.cell_range(cell, p, derivative);
                        result = Some(result.map_or(value, |old| hull(old, value)));
                    }
                }
            }
        }
        let result = result.unwrap_or(global);
        let lower = result.lo().max(global.lo());
        let upper = result.hi().min(global.hi());
        if lower <= upper { Interval::new(lower, upper) } else { global }
    }
}

impl CutSdf3 for SampledSdf3 {
    fn value(&self, p: [f64; 3]) -> f64 {
        if !self.covered(p, p) { return f64::NAN; }
        let cell: [usize; 3] = std::array::from_fn(|a|
            self.axes[a].partition_point(|&x| x <= p[a]).saturating_sub(1).min(self.shape[a] - 2));
        let mut values = self.corners(cell);
        for a in 0..3 {
            let left = self.axes[a][cell[a]];
            let t = ((p[a] - left) / (self.axes[a][cell[a] + 1] - left)).clamp(0.0, 1.0);
            let stride = 1 << a;
            for at in (0..8).step_by(2 * stride) {
                values[at] = values[at] * (1.0 - t) + values[at + stride] * t;
            }
        }
        values[0]
    }

    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval { self.range(lo, hi, None) }

    fn derivative_enclose(&self, lo: [f64; 3], hi: [f64; 3], axis: HeightAxis) -> Interval {
        let axis = match axis { HeightAxis::X => 0, HeightAxis::Y => 1, HeightAxis::Z => 2 };
        self.range(lo, hi, Some(axis))
    }
}

#[cfg(test)]
#[path = "sampled3_tests.rs"]
mod tests;

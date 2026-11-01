//! Ray-family clipping on the original cell hierarchy, with outward endpoints.
//!
//! For A in the start box and B in the end box, x(t)=A+t(B-A), 0<=t<=1.
//! Each strictly signed direction admits interval entry/exit parameters. A
//! zero-containing direction is never divided: its swept coordinate either
//! excludes the slab, is wholly inside it, or retains uncertain ownership.
//! In the last case the cell length has zero lower bound. Independent cell
//! bounds may overcount the upper integral, but enclose every member ray.

use super::*;

/// Version of the endpoint-box and uncertain-direction interpretation.
/// The original point-segment interpretation and arithmetic are unchanged.
pub const GRID_BUNDLE_ENCLOSURE_SEMANTICS_VERSION: u32 = 1;

impl AbsorbingGrid<'_> {
    /// Enclose column mass for all segments joining two finite coordinate boxes.
    /// Outside this original grid is vacuum. Degenerate member segments have
    /// zero column mass; a box pair containing both zero and nonzero lengths is
    /// not rejected or reduced to a nominal direction. Singleton boxes dispatch
    /// to the original endpoint calculation, preserving its bits and counters.
    /// The original cell-work limit and the additional node limit both apply.
    pub fn column_density_box_enclosure(
        &self, cx: &Cx<'_>, start: [Interval; 3], end: [Interval; 3],
        budget: GridEnclosureBudget,
    ) -> Result<ColumnDensityEnclosure, VolumeTransportError> {
        self.visit_box_column(cx, start, end, budget, |_, _| Ok(()))
    }

    pub(super) fn visit_box_column(
        &self, cx: &Cx<'_>, start: [Interval; 3], end: [Interval; 3],
        budget: GridEnclosureBudget,
        visit: impl FnMut(usize, Interval) -> Result<(), VolumeTransportError>,
    ) -> Result<ColumnDensityEnclosure, VolumeTransportError> {
        checkpoint(cx)?;
        if start.iter().chain(&end).any(|v| !finite(*v)) {
            return Err(VolumeTransportError::Ray);
        }
        if start.iter().chain(&end).all(|v| v.lo() == v.hi()) {
            let a = start.map(|v| v.lo());
            let b = end.map(|v| v.lo());
            return self.visit_column(cx, Point3::new(a[0], a[1], a[2]),
                Point3::new(b[0], b[1], b[2]), budget, visit);
        }
        let bundle = Bundle::new(start, end)?;
        self.visit_segments(cx, bundle.length, budget,
            |lower, upper| bundle.overlap(lower, upper), visit)
    }
}

struct Bundle {
    start: [Interval; 3],
    end: [Interval; 3],
    delta: [Interval; 3],
    length: Interval,
}
impl Bundle {
    fn new(start: [Interval; 3], end: [Interval; 3]) -> Result<Self, VolumeTransportError> {
        let delta = std::array::from_fn(|i| {
            if start[i].lo() == start[i].hi() && start[i] == end[i] {
                Interval::point(0.0)
            } else { end[i] - start[i] }
        });
        if delta.iter().any(|v| !finite(*v)) { return Err(VolumeTransportError::Unrepresentable); }
        let scale = delta.iter().map(|v| v.lo().abs().max(v.hi().abs())).fold(0.0_f64, f64::max);
        let length = if scale == 0.0 { Interval::point(0.0) } else {
            let square = delta.map(|v| {
                let x = v / Interval::point(scale);
                nonnegative(x * x)
            });
            nonnegative(nonnegative(square[0] + square[1] + square[2]).sqrt() * Interval::point(scale))
        };
        if !finite(length) { return Err(VolumeTransportError::Unrepresentable); }
        Ok(Self { start, end, delta, length })
    }

    fn overlap(&self, lower: [Interval; 3], upper: [Interval; 3]) -> Option<(Interval, bool)> {
        let mut enter = Interval::point(0.0);
        let mut exit = Interval::point(1.0);
        let mut uncertain = false;
        for axis in 0..3 {
            if self.delta[axis].contains_zero() {
                // Every coordinate along every member ray lies in this hull.
                let swept = self.start[axis].hull(self.end[axis]);
                if swept.hi() < lower[axis].lo() || swept.lo() >= upper[axis].hi() {
                    return None;
                }
                // A strict upper test preserves the native half-open ownership.
                uncertain |= swept.lo() < lower[axis].hi() || swept.hi() >= upper[axis].lo();
                continue;
            }
            let first = (lower[axis] - self.start[axis]) / self.delta[axis];
            let second = (upper[axis] - self.start[axis]) / self.delta[axis];
            let (near, far) = if self.delta[axis].lo() > 0.0 { (first, second) } else { (second, first) };
            enter = Interval::new(enter.lo().max(near.lo()), enter.hi().max(near.hi()));
            exit = Interval::new(exit.lo().min(far.lo()), exit.hi().min(far.hi()));
            if enter.lo() >= exit.hi() { return None; }
        }
        let width = nonnegative(exit - enter);
        let high = width.hi().min(1.0);
        if high == 0.0 { return None; }
        Some((Interval::new(if uncertain { 0.0 } else { width.lo() }, high), uncertain))
    }
}

#[cfg(test)]
mod tests;

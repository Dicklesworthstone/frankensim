//! Bounded evaluation of the native polar/azimuth interpolant, not a fixture fit.

use super::super::AngularTable;
use crate::photometry::PhotometryError;
use fs_exec::Cx;
use fs_ivl::{Interval, angle::atan2_upper_half_plane};
use std::f64::consts::TAU;

/// Maximum table cells examined for one source/probe numerical enclosure.
/// Larger direction uncertainty uses the valid global `[0, peak]` intensity
/// envelope instead of unbounded work or an unjustified midpoint evaluation.
/// Ordinary point directions visit only neighboring cells via binary searches.
pub const MAX_ANGULAR_CELLS_PER_ENCLOSURE: usize = 4096;

pub(super) fn table_enclosure(
    table: &AngularTable,
    peak: f64,
    local: [Interval; 3],
    mu: Interval,
    cx: &Cx<'_>,
) -> Result<Interval, PhotometryError> {
    cx.checkpoint()?;
    let global = Interval::new(0.0, peak);
    let arcs = if table.phi.len() == 2 {
        // The original two-column seam has equal values on each row. Preserve
        // its axisymmetric meaning without performing any inverse trig work.
        [Some(Interval::point(0.0)), None]
    } else {
        azimuth_arcs(local[0], local[1])
    };
    let first_row = table.cosine.partition_point(|v| *v > mu.hi())
        .saturating_sub(1).min(table.cosine.len() - 2);
    let last_row = table.cosine.partition_point(|v| *v >= mu.lo())
        .saturating_sub(1).min(table.cosine.len() - 2);
    let mut cells = 0usize;
    for arc in arcs.into_iter().flatten() {
        let (first, last) = columns(&table.phi, arc);
        cells = cells.saturating_add((last_row - first_row + 1).saturating_mul(last - first + 1));
    }
    if cells > MAX_ANGULAR_CELLS_PER_ENCLOSURE {
        return Ok(global);
    }
    let mut result: Option<Interval> = None;
    let mut visited = 0usize;
    for arc in arcs.into_iter().flatten() {
        let (first, last) = columns(&table.phi, arc);
        for row in first_row..=last_row {
            let lower = mu.lo().max(table.cosine[row + 1]);
            let upper = mu.hi().min(table.cosine[row]);
            if lower > upper { continue; }
            // t grows in the descending-cosine direction.
            let t = fraction(
                Interval::new(-upper, -lower),
                -table.cosine[row], -table.cosine[row + 1],
            );
            for column in first..=last {
                if visited % 64 == 0 { cx.checkpoint()?; }
                visited += 1;
                let lo = arc.lo().max(table.phi[column]);
                let hi = arc.hi().min(table.phi[column + 1]);
                if lo > hi { continue; }
                let u = fraction(Interval::new(lo, hi), table.phi[column], table.phi[column + 1]);
                let i = row * table.phi.len() + column;
                let a = Interval::point(table.values[i]);
                let b = Interval::point(table.values[i + 1]);
                let c = Interval::point(table.values[i + table.phi.len()]);
                let d = Interval::point(table.values[i + table.phi.len() + 1]);
                // A bilinear interpolant's extrema on a parameter rectangle
                // occur at its four corners, including saddle-shaped cells.
                for azimuth in [u.lo(), u.hi()] {
                    let top = convex(a, b, azimuth);
                    let bottom = convex(c, d, azimuth);
                    for polar in [t.lo(), t.hi()] {
                        let value = convex(top, bottom, polar);
                        result = Some(result.map_or(value, |old| old.hull(value)));
                    }
                }
            }
        }
    }
    cx.checkpoint()?;
    // All admitted knots lie in [0, peak]; convex interpolation preserves it.
    Ok(result.unwrap_or(global))
}

fn columns(phi: &[f64], arc: Interval) -> (usize, usize) {
    let first = phi.partition_point(|v| *v < arc.lo()).saturating_sub(1).min(phi.len() - 2);
    let last = phi.partition_point(|v| *v <= arc.hi()).saturating_sub(1).min(phi.len() - 2);
    (first, last)
}

/// Keep seam-adjacent directions as TWO arcs, never their full-circle hull.
/// Wrapping uses the native retained binary64 TAU; this encloses the declared
/// numerical source law rather than silently changing its azimuth coordinate.
fn azimuth_arcs(x: Interval, y: Interval) -> [Option<Interval>; 2] {
    if x.contains_zero() && y.contains_zero() {
        // Even a tiny box around a polar axis can include every azimuth. Pole
        // rows have identical values, so the table interpolation still tightens.
        return [Some(Interval::new(0.0, TAU)), None];
    }
    let mut result = [None, None];
    if y.hi() >= 0.0 {
        result[0] = atan2_upper_half_plane(Interval::new(y.lo().max(0.0), y.hi()), x)
            .map(clamp_azimuth);
    }
    if y.lo() < 0.0 {
        result[1] = atan2_upper_half_plane(Interval::new((-y.hi()).max(0.0), -y.lo()), x)
            .map(|angle| clamp_azimuth(Interval::point(TAU) - angle));
    }
    result
}

fn clamp_azimuth(v: Interval) -> Interval {
    Interval::new(v.lo().clamp(0.0, TAU), v.hi().clamp(0.0, TAU))
}

fn fraction(v: Interval, a: f64, b: f64) -> Interval {
    if v.lo() == a && v.hi() == a { return Interval::point(0.0); }
    if v.lo() == b && v.hi() == b { return Interval::point(1.0); }
    let denominator = Interval::point(b) - Interval::point(a);
    if denominator.lo() <= 0.0 { return Interval::new(0.0, 1.0); }
    let value = (v - Interval::point(a)) / denominator;
    Interval::new(value.lo().max(0.0), value.hi().min(1.0))
}

fn convex(a: Interval, b: Interval, t: f64) -> Interval {
    if t == 0.0 { return a; }
    if t == 1.0 { return b; }
    let t = Interval::point(t);
    let value = (Interval::point(1.0) - t) * a + t * b;
    // This also removes spurious finite overflows when both knots are near
    // f64::MAX, and preserves an exactly black cell instead of a tiny glow.
    Interval::new(value.lo().max(a.lo().min(b.lo())), value.hi().min(a.hi().max(b.hi())))
}

#[cfg(test)]
mod tests;

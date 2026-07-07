//! Cut-cell quadrature: accurate integration on implicitly-defined
//! subdomains, with error control by subdivision depth.
//!
//! The scheme is tessellation-based with CERTIFIED routing: a cut cell
//! recursively subdivides; every sub-cell is re-classified through the
//! SDF's interval enclosure, so fully-inside sub-cells get exact tensor
//! Gauss rules and fully-outside sub-cells vanish — only genuinely cut
//! sub-cells reach the marching-squares finale, where edge crossings
//! are located by bisection (resolution ~2⁻⁵⁰ of the sub-cell) and the
//! inside region is polygonized and integrated with a degree-2-exact
//! midpoint triangle rule.
//!
//! Error control: for a linear level set the crossings, the polygon,
//! and hence ALL quadratic moments are exact (the G0 exactness
//! battery). For curved interfaces the geometric error is
//! O((h/2^depth)²) per cell — quadratic convergence in depth,
//! measured and asserted by the conformance suite. Features smaller
//! than the finest subdivision that never change a corner sign can be
//! missed by the QUADRATURE (bounded by one sub-cell area) — never by
//! the CLASSIFICATION, which is interval-certified at every level; the
//! saddle/blob ambiguity spends a bounded extra-subdivision budget
//! before falling back to corner-sign integration.

use crate::sdf::CutSdf;

/// Bulk and interface rules for one cut cell, in global coordinates.
/// Bulk weights sum to the inside area (up to the documented error);
/// interface weights sum to the interface length, each point carrying
/// the OUTWARD unit normal (∇φ direction).
#[derive(Debug, Clone, Default)]
pub struct CutRules {
    /// Bulk points: (position, weight).
    pub bulk: Vec<([f64; 2], f64)>,
    /// Interface points: (position, weight, outward unit normal).
    pub iface: Vec<([f64; 2], f64, [f64; 2])>,
}

/// 3-point Gauss–Legendre on [-1, 1].
const G3: [(f64, f64); 3] = [
    (-0.774_596_669_241_483_4, 0.555_555_555_555_555_6),
    (0.0, 0.888_888_888_888_888_9),
    (0.774_596_669_241_483_4, 0.555_555_555_555_555_6),
];

/// Extra subdivision budget for saddle/blob ambiguity at depth 0.
const EXTRA: u32 = 2;

/// Push a 3×3 tensor Gauss rule for the full box (degree-5 exact per
/// axis — exact for every integrand this crate assembles).
pub fn tensor_gauss(lo: [f64; 2], hi: [f64; 2], out: &mut Vec<([f64; 2], f64)>) {
    let mx = 0.5 * (lo[0] + hi[0]);
    let my = 0.5 * (lo[1] + hi[1]);
    let sx = 0.5 * (hi[0] - lo[0]);
    let sy = 0.5 * (hi[1] - lo[1]);
    for &(gx, wx) in &G3 {
        for &(gy, wy) in &G3 {
            out.push(([mx + sx * gx, my + sy * gy], wx * wy * sx * sy));
        }
    }
}

/// Build the bulk + interface rules for one (certified-Cut) cell.
#[must_use]
pub fn cut_cell_rules(sdf: &dyn CutSdf, lo: [f64; 2], hi: [f64; 2], depth: u32) -> CutRules {
    let mut rules = CutRules::default();
    worker(sdf, lo, hi, depth, EXTRA, &mut rules);
    rules
}

fn worker(sdf: &dyn CutSdf, lo: [f64; 2], hi: [f64; 2], depth: u32, extra: u32, out: &mut CutRules) {
    let iv = sdf.enclose(lo, hi);
    if iv.hi() < 0.0 {
        tensor_gauss(lo, hi, &mut out.bulk);
        return;
    }
    if iv.lo() > 0.0 {
        return;
    }
    if depth > 0 {
        recurse(sdf, lo, hi, depth - 1, extra, out);
        return;
    }
    // Finest level: marching squares on the corner signs.
    let corners = [
        [lo[0], lo[1]],
        [hi[0], lo[1]],
        [hi[0], hi[1]],
        [lo[0], hi[1]],
    ];
    let phi: Vec<f64> = corners.iter().map(|&p| sdf.value(p)).collect();
    let inside: Vec<bool> = phi.iter().map(|&v| v <= 0.0).collect();
    let mut crossings: [Option<[f64; 2]>; 4] = [None; 4];
    let mut ncross = 0;
    for e in 0..4 {
        if inside[e] != inside[(e + 1) % 4] {
            crossings[e] = Some(bisect_crossing(sdf, corners[e], corners[(e + 1) % 4]));
            ncross += 1;
        }
    }
    if ncross == 0 {
        // Interval says "maybe cut" but no corner sign change: a
        // sub-resolution feature (tangency or interior blob). Spend the
        // extra budget, then fall back to the corner-sign picture.
        if extra > 0 {
            recurse(sdf, lo, hi, 0, extra - 1, out);
        } else if inside[0] {
            tensor_gauss(lo, hi, &mut out.bulk);
        }
        return;
    }
    if ncross == 4 {
        // Saddle: prefer more resolution; if exhausted, resolve the
        // connectivity by the center sign.
        if extra > 0 {
            recurse(sdf, lo, hi, 0, extra - 1, out);
            return;
        }
        saddle_rules(sdf, &corners, &inside, &crossings, out);
        return;
    }
    // Regular case (2 crossings): walk the boundary counterclockwise,
    // emitting inside corners and crossings — a simple polygon.
    let mut poly: Vec<[f64; 2]> = Vec::with_capacity(6);
    for e in 0..4 {
        if inside[e] {
            poly.push(corners[e]);
        }
        if let Some(x) = crossings[e] {
            poly.push(x);
        }
    }
    polygon_rule(&poly, &mut out.bulk);
    let xs: Vec<[f64; 2]> = crossings.iter().flatten().copied().collect();
    chord_rule(sdf, xs[0], xs[1], &mut out.iface);
}

fn recurse(sdf: &dyn CutSdf, lo: [f64; 2], hi: [f64; 2], depth: u32, extra: u32, out: &mut CutRules) {
    let mx = 0.5 * (lo[0] + hi[0]);
    let my = 0.5 * (lo[1] + hi[1]);
    worker(sdf, lo, [mx, my], depth, extra, out);
    worker(sdf, [mx, lo[1]], [hi[0], my], depth, extra, out);
    worker(sdf, [lo[0], my], [mx, hi[1]], depth, extra, out);
    worker(sdf, [mx, my], hi, depth, extra, out);
}

/// Saddle finale: the center sign decides which crossing pairs bound
/// the inside region. Center inside → the hexagon walk (both chords
/// cut off outside corner lobes); center outside → two inside corner
/// triangles.
fn saddle_rules(
    sdf: &dyn CutSdf,
    corners: &[[f64; 2]; 4],
    inside: &[bool],
    crossings: &[Option<[f64; 2]>; 4],
    out: &mut CutRules,
) {
    let center = [
        0.5 * (corners[0][0] + corners[2][0]),
        0.5 * (corners[0][1] + corners[2][1]),
    ];
    let center_in = sdf.value(center) <= 0.0;
    if center_in {
        let mut poly: Vec<[f64; 2]> = Vec::with_capacity(6);
        for e in 0..4 {
            if inside[e] {
                poly.push(corners[e]);
            }
            if let Some(x) = crossings[e] {
                poly.push(x);
            }
        }
        polygon_rule(&poly, &mut out.bulk);
        for k in 0..4 {
            if !inside[k] {
                let a = crossings[(k + 3) % 4].expect("saddle has all crossings");
                let b = crossings[k].expect("saddle has all crossings");
                chord_rule(sdf, a, b, &mut out.iface);
            }
        }
    } else {
        for k in 0..4 {
            if inside[k] {
                let a = crossings[(k + 3) % 4].expect("saddle has all crossings");
                let b = crossings[k].expect("saddle has all crossings");
                polygon_rule(&[a, corners[k], b], &mut out.bulk);
                chord_rule(sdf, a, b, &mut out.iface);
            }
        }
    }
}

/// Bisection for the interface crossing on a segment with a corner
/// sign change (~2⁻⁵⁰ of the segment; exact for linear φ up to
/// roundoff).
fn bisect_crossing(sdf: &dyn CutSdf, a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    let sa = sdf.value(a) > 0.0;
    let (mut t0, mut t1) = (0.0f64, 1.0f64);
    for _ in 0..50 {
        let tm = 0.5 * (t0 + t1);
        let p = [a[0] + tm * (b[0] - a[0]), a[1] + tm * (b[1] - a[1])];
        if (sdf.value(p) > 0.0) == sa {
            t0 = tm;
        } else {
            t1 = tm;
        }
    }
    let t = 0.5 * (t0 + t1);
    [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]
}

/// Fan-triangulate a simple polygon (walk order) and push the
/// degree-2-exact midpoint rule per triangle. Signed areas make the
/// fan correct for non-convex walk polygons.
fn polygon_rule(poly: &[[f64; 2]], out: &mut Vec<([f64; 2], f64)>) {
    if poly.len() < 3 {
        return;
    }
    for k in 1..poly.len() - 1 {
        let (p, q, r) = (poly[0], poly[k], poly[k + 1]);
        let signed_area =
            0.5 * ((q[0] - p[0]) * (r[1] - p[1]) - (r[0] - p[0]) * (q[1] - p[1]));
        if signed_area == 0.0 {
            continue;
        }
        let w = signed_area / 3.0;
        out.push(([0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1])], w));
        out.push(([0.5 * (q[0] + r[0]), 0.5 * (q[1] + r[1])], w));
        out.push(([0.5 * (r[0] + p[0]), 0.5 * (r[1] + p[1])], w));
    }
}

/// 2-point Gauss along an interface chord; weights carry the chord
/// length, normals come from the (normalized) SDF gradient — outward
/// of Ω by the negative-inside convention.
fn chord_rule(sdf: &dyn CutSdf, a: [f64; 2], b: [f64; 2], out: &mut Vec<([f64; 2], f64, [f64; 2])>) {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let len = (dx * dx + dy * dy).sqrt();
    if len == 0.0 {
        return;
    }
    let g = 0.5 / 3.0f64.sqrt();
    for t in [0.5 - g, 0.5 + g] {
        let p = [a[0] + t * dx, a[1] + t * dy];
        let grad = sdf.gradient(p);
        let gn = (grad[0] * grad[0] + grad[1] * grad[1]).sqrt();
        let normal = if gn > 1e-300 {
            [grad[0] / gn, grad[1] / gn]
        } else {
            // Degenerate gradient: fall back to the chord perpendicular,
            // signed by a φ probe.
            let perp = [dy / len, -dx / len];
            let eps = 1e-9 * len.max(1e-12);
            let plus = sdf.value([p[0] + eps * perp[0], p[1] + eps * perp[1]]);
            let minus = sdf.value([p[0] - eps * perp[0], p[1] - eps * perp[1]]);
            if plus >= minus {
                perp
            } else {
                [-perp[0], -perp[1]]
            }
        };
        out.push((p, 0.5 * len, normal));
    }
}

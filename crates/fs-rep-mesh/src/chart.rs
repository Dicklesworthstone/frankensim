//! The mesh chart: signed distance = exact point-triangle distance (BVH-
//! accelerated) with the generalized-winding-number SIGN — robust even on
//! imperfect input — plus watertight ray-triangle intersection (Woop et
//! al. 2013, the convention LUMEN's tracer will share).
//!
//! Authority honesty: unsigned distance to a triangle set is 1-Lipschitz, but
//! multiplying it by a generalized-winding sign on an arbitrary open,
//! self-intersecting, or inconsistently oriented soup need not preserve that
//! theorem or define the intended solid. The raw chart therefore exposes
//! nominal `Estimate` evidence and no certified Lipschitz bound. A future
//! validated wrapper may upgrade both only after sufficient solid-validity,
//! orientation, nesting, and self-intersection proofs.

use crate::winding::{Soup, WindingOctree};
use fs_evidence::NumericalCertificate;
use fs_exec::Cx;
use fs_geom::{Aabb, Chart, ChartSample, Differentiability, Point3, Vec3};

/// Exact closest distance from `p` to segment `[s0, s1]`. A zero-length
/// segment (`s0 == s1`) degrades gracefully to the point distance (`t = 0`).
fn point_segment_distance(p: Point3, s0: Point3, s1: Point3) -> f64 {
    let e = s1.delta_from(s0);
    let ee = e.dot(e);
    let t = if ee > 0.0 {
        (p.delta_from(s0).dot(e) / ee).clamp(0.0, 1.0)
    } else {
        0.0
    };
    p.delta_from(s0.offset(e.scale(t))).norm()
}

/// Exact closest distance from `p` to triangle `[a, b, c]` (Ericson,
/// Real-Time Collision Detection §5.1.5 — the standard region test).
#[must_use]
pub fn point_triangle_distance(p: Point3, a: Point3, b: Point3, c: Point3) -> f64 {
    let ab = b.delta_from(a);
    let ac = c.delta_from(a);
    // Ericson's region test divides by |ab|², |ac|², |bc|², and (twice) the
    // triangle area — every one of which is zero EXACTLY when the triangle is
    // degenerate (a repeated vertex, or three collinear vertices). Their shared
    // witness is |ab×ac|² = |ab|²·|ac|² − (ab·ac)² (Lagrange): strictly > 0 for a
    // proper triangle, ≤ 0 (zero, modulo cancellation on a razor-thin one) when
    // degenerate. In that case the "triangle" collapses to a segment/point, so
    // the true closest distance is the nearest of the three edge distances —
    // returning it keeps the result well-defined (never a 0/0 NaN, e.g. the
    // edge-AB branch's `d1/(d1-d3)` with `a == b`), as the CONTRACT promises.
    // Non-degenerate triangles skip this and hit the unchanged region test.
    let dab_ac = ab.dot(ac);
    if ab.dot(ab) * ac.dot(ac) - dab_ac * dab_ac <= 0.0 {
        return point_segment_distance(p, a, b)
            .min(point_segment_distance(p, b, c))
            .min(point_segment_distance(p, a, c));
    }
    let ap = p.delta_from(a);
    let d1 = ab.dot(ap);
    let d2 = ac.dot(ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return ap.norm();
    }
    let bp = p.delta_from(b);
    let d3 = ab.dot(bp);
    let d4 = ac.dot(bp);
    if d3 >= 0.0 && d4 <= d3 {
        return bp.norm();
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return p.delta_from(a.offset(ab.scale(v))).norm();
    }
    let cp = p.delta_from(c);
    let d5 = ab.dot(cp);
    let d6 = ac.dot(cp);
    if d6 >= 0.0 && d5 <= d6 {
        return cp.norm();
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return p.delta_from(a.offset(ac.scale(w))).norm();
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        let bc = c.delta_from(b);
        return p.delta_from(b.offset(bc.scale(w))).norm();
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    let closest = a.offset(ab.scale(v)).offset(ac.scale(w));
    p.delta_from(closest).norm()
}

/// Watertight ray-triangle intersection (Woop, Benthin, Wald 2013):
/// dimension permutation + shear puts the ray on +z, then signed edge
/// functions with consistent orientation guarantee shared edges/vertices
/// never leak or double-hit. Returns the hit parameter `t ≥ 0`.
#[must_use]
pub fn ray_triangle_watertight(
    origin: Point3,
    dir: Vec3,
    a: Point3,
    b: Point3,
    c: Point3,
) -> Option<f64> {
    // Permute so |dir.kz| is the largest component.
    let abs = [dir.x.abs(), dir.y.abs(), dir.z.abs()];
    let kz = if abs[0] > abs[1] && abs[0] > abs[2] {
        0
    } else if abs[1] > abs[2] {
        1
    } else {
        2
    };
    let kx = (kz + 1) % 3;
    let ky = (kz + 2) % 3;
    let comp = |v: Vec3, k: usize| match k {
        0 => v.x,
        1 => v.y,
        _ => v.z,
    };
    let dz = comp(dir, kz);
    let (kx, ky) = if dz < 0.0 { (ky, kx) } else { (kx, ky) };
    // Shear constants.
    let sx = comp(dir, kx) / dz;
    let sy = comp(dir, ky) / dz;
    let sz = 1.0 / dz;
    // Vertices relative to the origin, sheared.
    let rel = |v: Point3| v.delta_from(origin);
    let (va, vb, vc) = (rel(a), rel(b), rel(c));
    let ax = comp(va, kx) - sx * comp(va, kz);
    let ay = comp(va, ky) - sy * comp(va, kz);
    let bx = comp(vb, kx) - sx * comp(vb, kz);
    let by = comp(vb, ky) - sy * comp(vb, kz);
    let cx = comp(vc, kx) - sx * comp(vc, kz);
    let cy = comp(vc, ky) - sy * comp(vc, kz);
    // Signed edge functions.
    let u = cx * by - cy * bx;
    let v = ax * cy - ay * cx;
    let w = bx * ay - by * ax;
    // Fall back to f64-widened re-evaluation on exact-zero edges (we are
    // already in f64; the paper's double-fallback becomes an exact-tie
    // rule here).
    if (u < 0.0 || v < 0.0 || w < 0.0) && (u > 0.0 || v > 0.0 || w > 0.0) {
        return None;
    }
    let det = u + v + w;
    if det == 0.0 {
        return None;
    }
    // Sheared z for the hit distance.
    let az = sz * comp(va, kz);
    let bz = sz * comp(vb, kz);
    let cz = sz * comp(vc, kz);
    // Backface-agnostic forward hits only. Shared edges/vertices never
    // LEAK (zero edge functions with consistent signs still hit); an
    // exact edge ray may double-count across the two incident triangles —
    // fine for nearest-hit queries, documented for parity users.
    let t = (u * az + v * bz + w * cz) / det;
    if t >= 0.0 { Some(t) } else { None }
}

/// A flat median-split AABB BVH over soup triangles (shared design target
/// with LUMEN's tracer — one implementation, two consumers).
pub struct Bvh {
    nodes: Vec<BvhNode>,
    /// Triangle order after building (leaves index into this).
    order: Vec<u32>,
}

struct BvhNode {
    bounds: Aabb,
    /// Leaf: (start, count) into `order`; internal: (left child, right
    /// child) with count == 0 marker.
    a: u32,
    b: u32,
    leaf_count: u32,
}

impl Bvh {
    /// Build over a soup.
    #[must_use]
    pub fn build(soup: &Soup) -> Self {
        let mut order: Vec<u32> = (0..soup.triangles.len() as u32).collect();
        let mut nodes = Vec::new();
        let n = order.len();
        Self::build_range(soup, &mut order, 0, n, &mut nodes);
        Bvh { nodes, order }
    }

    fn tri_bounds(soup: &Soup, t: u32) -> Aabb {
        let [a, b, c] = soup.tri(t as usize);
        Aabb::new(a, b).union(&Aabb::new(c, c))
    }

    fn centroid(soup: &Soup, t: u32) -> Point3 {
        let [a, b, c] = soup.tri(t as usize);
        Point3::new(
            (a.x + b.x + c.x) / 3.0,
            (a.y + b.y + c.y) / 3.0,
            (a.z + b.z + c.z) / 3.0,
        )
    }

    fn build_range(
        soup: &Soup,
        order: &mut [u32],
        start: usize,
        count: usize,
        nodes: &mut Vec<BvhNode>,
    ) -> u32 {
        let mut bounds = Self::tri_bounds(soup, order[start]);
        for &t in &order[start..start + count] {
            bounds = bounds.union(&Self::tri_bounds(soup, t));
        }
        let idx = nodes.len() as u32;
        if count <= 4 {
            nodes.push(BvhNode {
                bounds,
                a: start as u32,
                b: 0,
                leaf_count: count as u32,
            });
            return idx;
        }
        // Median split along the widest axis of the centroid bounds.
        let extent = [
            bounds.max.x - bounds.min.x,
            bounds.max.y - bounds.min.y,
            bounds.max.z - bounds.min.z,
        ];
        let axis = if extent[0] >= extent[1] && extent[0] >= extent[2] {
            0
        } else if extent[1] >= extent[2] {
            1
        } else {
            2
        };
        let slice = &mut order[start..start + count];
        slice.sort_unstable_by(|&x, &y| {
            let cx = Self::centroid(soup, x);
            let cy = Self::centroid(soup, y);
            let (vx, vy) = match axis {
                0 => (cx.x, cy.x),
                1 => (cx.y, cy.y),
                _ => (cx.z, cy.z),
            };
            vx.total_cmp(&vy).then(x.cmp(&y)) // deterministic ties
        });
        let mid = count / 2;
        nodes.push(BvhNode {
            bounds,
            a: 0,
            b: 0,
            leaf_count: 0,
        });
        let left = Self::build_range(soup, order, start, mid, nodes);
        let right = Self::build_range(soup, order, start + mid, count - mid, nodes);
        nodes[idx as usize].a = left;
        nodes[idx as usize].b = right;
        idx
    }

    fn aabb_dist_sq(bounds: &Aabb, p: Point3) -> f64 {
        let dx = (bounds.min.x - p.x).max(0.0).max(p.x - bounds.max.x);
        let dy = (bounds.min.y - p.y).max(0.0).max(p.y - bounds.max.y);
        let dz = (bounds.min.z - p.z).max(0.0).max(p.z - bounds.max.z);
        dx * dx + dy * dy + dz * dz
    }

    /// Closest distance from `p` to any triangle (branch-and-bound).
    #[must_use]
    pub fn closest_distance(&self, soup: &Soup, p: Point3) -> f64 {
        let mut best = f64::INFINITY;
        let mut stack = vec![0u32];
        while let Some(ni) = stack.pop() {
            let node = &self.nodes[ni as usize];
            if Self::aabb_dist_sq(&node.bounds, p) >= best * best {
                continue;
            }
            if node.leaf_count > 0 {
                for i in node.a..node.a + node.leaf_count {
                    let t = self.order[i as usize];
                    let [a, b, c] = soup.tri(t as usize);
                    best = best.min(point_triangle_distance(p, a, b, c));
                }
            } else {
                stack.push(node.a);
                stack.push(node.b);
            }
        }
        best
    }

    /// Nearest forward ray hit over the soup (watertight per triangle).
    #[must_use]
    pub fn raycast(&self, soup: &Soup, origin: Point3, dir: Vec3, t_max: f64) -> Option<f64> {
        let mut best: Option<f64> = None;
        let mut stack = vec![0u32];
        while let Some(ni) = stack.pop() {
            let node = &self.nodes[ni as usize];
            // Conservative slab test.
            if !ray_hits_aabb(&node.bounds, origin, dir, best.unwrap_or(t_max)) {
                continue;
            }
            if node.leaf_count > 0 {
                for i in node.a..node.a + node.leaf_count {
                    let t = self.order[i as usize];
                    let [a, b, c] = soup.tri(t as usize);
                    if let Some(hit) = ray_triangle_watertight(origin, dir, a, b, c)
                        && hit <= best.unwrap_or(t_max)
                    {
                        best = Some(hit);
                    }
                }
            } else {
                stack.push(node.a);
                stack.push(node.b);
            }
        }
        best
    }
}

fn ray_hits_aabb(bounds: &Aabb, origin: Point3, dir: Vec3, t_max: f64) -> bool {
    let mut t0 = 0.0f64;
    let mut t1 = t_max;
    for axis in 0..3 {
        let (o, d, lo, hi) = match axis {
            0 => (origin.x, dir.x, bounds.min.x, bounds.max.x),
            1 => (origin.y, dir.y, bounds.min.y, bounds.max.y),
            _ => (origin.z, dir.z, bounds.min.z, bounds.max.z),
        };
        if d.abs() < 1e-300 {
            if o < lo || o > hi {
                return false;
            }
            continue;
        }
        let inv = 1.0 / d;
        let (mut ta, mut tb) = ((lo - o) * inv, (hi - o) * inv);
        if ta > tb {
            core::mem::swap(&mut ta, &mut tb);
        }
        t0 = t0.max(ta);
        t1 = t1.min(tb);
        if t0 > t1 {
            return false;
        }
    }
    true
}

/// The soup-backed chart: robust signed distance for possibly imperfect
/// meshes (plan §7.2's mesh chart).
pub struct MeshChart {
    soup: Soup,
    bvh: Bvh,
    winding: WindingOctree,
    support: Aabb,
}

impl MeshChart {
    /// Build from a soup (BVH + winding octree constructed ONCE).
    #[must_use]
    pub fn new(soup: Soup) -> Self {
        let bvh = Bvh::build(&soup);
        let winding = WindingOctree::build(&soup, 2.0);
        let mut support = Aabb::new(soup.positions[0], soup.positions[0]);
        for &p in &soup.positions {
            support = support.union(&Aabb::new(p, p));
        }
        MeshChart {
            soup,
            bvh,
            winding,
            support,
        }
    }

    /// The triangle soup.
    #[must_use]
    pub fn soup(&self) -> &Soup {
        &self.soup
    }

    /// Watertight raycast against the mesh.
    #[must_use]
    pub fn raycast(&self, origin: Point3, dir: Vec3, t_max: f64) -> Option<f64> {
        self.bvh.raycast(&self.soup, origin, dir, t_max)
    }
}

impl Chart for MeshChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let dist = self.bvh.closest_distance(&self.soup, x);
        let sign = if self.winding.winding(&self.soup, x) > 0.5 {
            -1.0
        } else {
            1.0
        };
        let sd = sign * dist;
        // This is only a nominal accuracy scale. The unsigned magnitude is a
        // distance-to-set, but the generalized-winding sign can jump on an
        // invalid soup, so neither global Lipschitz nor abstract-region
        // distance authority follows from the raw representation.
        let eps = 1e-12 * (1.0 + x.delta_from(self.support.min).norm());
        let error = if sd.is_finite() && eps.is_finite() {
            NumericalCertificate::estimate(sd - eps, sd + eps)
        } else {
            NumericalCertificate::no_claim()
        };
        ChartSample {
            signed_distance: sd,
            gradient: None,
            lipschitz: None,
            error,
        }
    }

    fn support(&self) -> Aabb {
        self.support
    }

    fn name(&self) -> &'static str {
        "rep-mesh/soup"
    }

    fn differentiability(&self) -> Differentiability {
        Differentiability::C0
    }
}

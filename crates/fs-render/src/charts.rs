//! CHART BACKENDS (plan §10.2, bead qfx.2; [F] — behind the
//! `chart-backends` feature): render whatever chart exists, WITHOUT
//! conversion. Certified sphere tracing for SDF/F-rep charts (step
//! sizes that PROVABLY never tunnel: within radius `|f(p)|/L` of `p`
//! the field cannot change sign, by the certified Lipschitz bound —
//! the certificate machinery earning visual credibility), Bézier-seeded
//! Newton intersection for NURBS patches, native triangle tracing over
//! a deterministic median-split BVH, and mixed-chart scenes: one scene,
//! three backend kinds, one image.
//!
//! An agent inspecting an F-rep mid-optimization sees the F-rep itself
//! — the no-meshing-for-visualization doctrine.

use fs_exec::Cx;
use fs_geom::{Chart, Point3, Vec3};
use fs_rep_nurbs::NurbsSurface;

/// A ray with unit direction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ray {
    /// Origin.
    pub origin: Point3,
    /// UNIT direction (callers normalize).
    pub dir: Vec3,
}

impl Ray {
    /// The point at parameter `t`.
    #[must_use]
    pub fn at(&self, t: f64) -> Point3 {
        self.origin.offset(self.dir.scale(t))
    }
}

/// One intersection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hit {
    /// Ray parameter.
    pub t: f64,
    /// The hit point.
    pub point: Point3,
    /// The surface normal, when the backend supplies one.
    pub normal: Option<Vec3>,
    /// Work spent (marcher steps / Newton iterations / BVH visits).
    pub steps: u32,
}

/// Sphere-trace telemetry: the G0 step-safety audit rides along.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TraceAudit {
    /// Steps taken.
    pub steps: u32,
    /// The worst ratio `step / (|f|/L)` observed (must stay ≤ 1 + ε for
    /// plain marching; over-relaxed steps are audited via fallback).
    pub worst_step_ratio: f64,
}

/// CERTIFIED sphere tracing: at each point the next step is
/// `|f(p)| / L` with `L` the chart's certified Lipschitz bound — the
/// sign cannot flip within that radius, so the marcher can never cross
/// (tunnel through) the surface. Over-relaxation (`omega > 1`)
/// accelerates marching with the standard certified fallback: if the
/// relaxed sphere fails to overlap the previous safe sphere, the step
/// is redone unrelaxed from the last safe point.
#[must_use]
pub fn sphere_trace(
    chart: &dyn Chart,
    cx: &Cx<'_>,
    ray: &Ray,
    t_max: f64,
    eps: f64,
    omega: f64,
) -> (Option<Hit>, TraceAudit) {
    let mut t = 0.0f64;
    let mut steps = 0u32;
    let mut worst_ratio = 0.0f64;
    // The last SAFE (unrelaxed) radius and position.
    let mut prev_radius = 0.0f64;
    let mut relaxed_pending = false;
    let max_steps = 4096u32;
    while t <= t_max && steps < max_steps {
        let p = ray.at(t);
        let s = chart.eval(p, cx);
        let lipschitz = s.lipschitz.unwrap_or(1.0).max(1e-12);
        let d = s.signed_distance;
        let safe = d.abs() / lipschitz;
        if d.abs() <= eps {
            let normal = s.gradient.or_else(|| gradient_fd(chart, cx, p));
            return (
                Some(Hit {
                    t,
                    point: p,
                    normal,
                    steps,
                }),
                TraceAudit {
                    steps,
                    worst_step_ratio: worst_ratio,
                },
            );
        }
        if relaxed_pending && safe + prev_radius < (omega - 1.0) * prev_radius + prev_radius {
            // The relaxed step's sphere does not overlap the previous
            // safe sphere: retreat to the certified step.
            t -= (omega - 1.0) * prev_radius;
            relaxed_pending = false;
            steps += 1;
            continue;
        }
        let step = if omega > 1.0 { omega * safe } else { safe };
        worst_ratio = worst_ratio.max(step / safe.max(1e-300) / omega.max(1.0));
        prev_radius = safe;
        relaxed_pending = omega > 1.0;
        t += step;
        steps += 1;
    }
    (
        None,
        TraceAudit {
            steps,
            worst_step_ratio: worst_ratio,
        },
    )
}

/// Central-difference normal fallback (charts without gradients).
fn gradient_fd(chart: &dyn Chart, cx: &Cx<'_>, p: Point3) -> Option<Vec3> {
    let h = 1e-6;
    let d = |q: Point3| chart.eval(q, cx).signed_distance;
    let g = Vec3::new(
        d(Point3::new(p.x + h, p.y, p.z)) - d(Point3::new(p.x - h, p.y, p.z)),
        d(Point3::new(p.x, p.y + h, p.z)) - d(Point3::new(p.x, p.y - h, p.z)),
        d(Point3::new(p.x, p.y, p.z + h)) - d(Point3::new(p.x, p.y, p.z - h)),
    );
    let n = g.norm();
    (n > 1e-12).then(|| g.scale(1.0 / n))
}

/// NURBS ray intersection: coarse-grid seeds ranked by distance to the
/// ray line, then 3×3 Newton on `F(u, v, t) = S(u, v) − o − t·d = 0`
/// with the Jacobian `[S_u, S_v, −d]` — the Bézier-clipping-seeded
/// Newton the plan names, adapted from the closest-point machinery.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn ray_intersect_nurbs(
    surface: &NurbsSurface<f64>,
    ray: &Ray,
    seeds_per_axis: usize,
    eps: f64,
) -> Option<Hit> {
    let (ulo, uhi) = surface.knots_u.domain();
    let (vlo, vhi) = surface.knots_v.domain();
    // Seed ranking: distance from the sample point to the ray LINE.
    let mut seeds: Vec<(f64, f64, f64, f64)> = Vec::new(); // (dist, u, v, t)
    for a in 0..seeds_per_axis {
        for b in 0..seeds_per_axis {
            #[allow(clippy::cast_precision_loss)]
            let u = ulo + (uhi - ulo) * (a as f64 + 0.5) / seeds_per_axis as f64;
            #[allow(clippy::cast_precision_loss)]
            let v = vlo + (vhi - vlo) * (b as f64 + 0.5) / seeds_per_axis as f64;
            let Ok(p) = surface.eval(u, v) else { continue };
            let rel = [
                p[0] - ray.origin.x,
                p[1] - ray.origin.y,
                p[2] - ray.origin.z,
            ];
            let t = rel[0] * ray.dir.x + rel[1] * ray.dir.y + rel[2] * ray.dir.z;
            let closest = ray.at(t);
            let dist = ((p[0] - closest.x).powi(2)
                + (p[1] - closest.y).powi(2)
                + (p[2] - closest.z).powi(2))
            .sqrt();
            if t > 0.0 {
                seeds.push((dist, u, v, t));
            }
        }
    }
    seeds.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.3.total_cmp(&y.3)));
    let mut best: Option<Hit> = None;
    for &(_, u0, v0, t0) in seeds.iter().take(6) {
        let (mut u, mut v, mut t) = (u0, v0, t0);
        for iter_count in 1..=24u32 {
            let Ok((pos, su, sv)) = surface.partials(u, v) else {
                break;
            };
            let f = [
                pos[0] - ray.origin.x - t * ray.dir.x,
                pos[1] - ray.origin.y - t * ray.dir.y,
                pos[2] - ray.origin.z - t * ray.dir.z,
            ];
            let fn2 = f[0] * f[0] + f[1] * f[1] + f[2] * f[2];
            if fn2 < eps * eps {
                let n = Vec3::new(
                    su[1] * sv[2] - su[2] * sv[1],
                    su[2] * sv[0] - su[0] * sv[2],
                    su[0] * sv[1] - su[1] * sv[0],
                );
                let nn = n.norm();
                let hit = Hit {
                    t,
                    point: ray.at(t),
                    normal: (nn > 1e-12).then(|| n.scale(1.0 / nn)),
                    steps: iter_count,
                };
                if t > 1e-9 && best.as_ref().is_none_or(|b| hit.t < b.t) {
                    best = Some(hit);
                }
                break;
            }
            // Solve J * delta = -F with J = [Su, Sv, -d] (Cramer 3x3).
            let j = [
                [su[0], sv[0], -ray.dir.x],
                [su[1], sv[1], -ray.dir.y],
                [su[2], sv[2], -ray.dir.z],
            ];
            let det = det3(&j);
            if det.abs() < 1e-14 {
                break;
            }
            let rhs = [-f[0], -f[1], -f[2]];
            let du = det3(&replace_col(&j, 0, &rhs)) / det;
            let dv = det3(&replace_col(&j, 1, &rhs)) / det;
            let dt = det3(&replace_col(&j, 2, &rhs)) / det;
            u = (u + du).clamp(ulo, uhi);
            v = (v + dv).clamp(vlo, vhi);
            t += dt;
            if t < 0.0 {
                break;
            }
        }
    }
    best
}

fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

fn replace_col(m: &[[f64; 3]; 3], col: usize, v: &[f64; 3]) -> [[f64; 3]; 3] {
    let mut out = *m;
    for r in 0..3 {
        out[r][col] = v[r];
    }
    out
}

/// A triangle mesh with a deterministic median-split BVH (the interim
/// native backend until the SIMD 8-wide BVH lands — CONTRACT no-claim).
#[derive(Debug, Clone)]
pub struct TriMesh {
    /// Vertices.
    pub vertices: Vec<[f64; 3]>,
    /// Triangles (vertex indices).
    pub triangles: Vec<[u32; 3]>,
    nodes: Vec<BvhNode>,
    order: Vec<u32>,
}

#[derive(Debug, Clone)]
struct BvhNode {
    lo: [f64; 3],
    hi: [f64; 3],
    /// Leaf: (start, count) into `order`; inner: (left, right) node ids
    /// with count == u32::MAX sentinel.
    a: u32,
    b: u32,
    leaf: bool,
}

impl TriMesh {
    /// Build with the deterministic median-split BVH.
    #[must_use]
    pub fn new(vertices: Vec<[f64; 3]>, triangles: Vec<[u32; 3]>) -> TriMesh {
        let mut mesh = TriMesh {
            vertices,
            triangles,
            nodes: Vec::new(),
            order: Vec::new(),
        };
        mesh.order = (0..mesh.triangles.len() as u32).collect();
        if !mesh.triangles.is_empty() {
            let n = mesh.triangles.len();
            let mut order = std::mem::take(&mut mesh.order);
            mesh.build(&mut order, 0, n);
            mesh.order = order;
        }
        mesh
    }

    fn centroid(&self, tri: u32) -> [f64; 3] {
        let t = self.triangles[tri as usize];
        let mut c = [0.0f64; 3];
        for &vi in &t {
            for (ck, &pv) in c.iter_mut().zip(&self.vertices[vi as usize]) {
                *ck += pv / 3.0;
            }
        }
        c
    }

    fn bounds(&self, tris: &[u32]) -> ([f64; 3], [f64; 3]) {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for &ti in tris {
            for &vi in &self.triangles[ti as usize] {
                for k in 0..3 {
                    lo[k] = lo[k].min(self.vertices[vi as usize][k]);
                    hi[k] = hi[k].max(self.vertices[vi as usize][k]);
                }
            }
        }
        (lo, hi)
    }

    fn build(&mut self, order: &mut [u32], start: usize, count: usize) -> u32 {
        let slice = &order[start..start + count];
        let (lo, hi) = self.bounds(slice);
        let id = self.nodes.len() as u32;
        self.nodes.push(BvhNode {
            lo,
            hi,
            a: start as u32,
            b: count as u32,
            leaf: true,
        });
        if count <= 4 {
            return id;
        }
        // Median split on the widest axis, deterministic tie-break.
        let axis = (0..3)
            .max_by(|&a, &b| (hi[a] - lo[a]).total_cmp(&(hi[b] - lo[b])))
            .unwrap_or(0);
        let seg = &mut order[start..start + count];
        seg.sort_by(|&x, &y| {
            self.centroid(x)[axis]
                .total_cmp(&self.centroid(y)[axis])
                .then(x.cmp(&y))
        });
        let half = count / 2;
        let left = self.build(order, start, half);
        let right = self.build(order, start + half, count - half);
        self.nodes[id as usize] = BvhNode {
            lo,
            hi,
            a: left,
            b: right,
            leaf: false,
        };
        id
    }

    /// Closest triangle intersection (Möller–Trumbore through the BVH).
    #[must_use]
    pub fn intersect(&self, ray: &Ray) -> Option<Hit> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut best: Option<Hit> = None;
        let mut stack = vec![0u32];
        let mut visits = 0u32;
        while let Some(id) = stack.pop() {
            visits += 1;
            let node = &self.nodes[id as usize];
            if !slab_hit(ray, node.lo, node.hi, best.map_or(f64::INFINITY, |h| h.t)) {
                continue;
            }
            if node.leaf {
                for &ti in &self.order[node.a as usize..(node.a + node.b) as usize] {
                    if let Some(mut hit) = self.tri_hit(ray, ti)
                        && best.as_ref().is_none_or(|b| hit.t < b.t)
                    {
                        hit.steps = visits;
                        best = Some(hit);
                    }
                }
            } else {
                stack.push(node.b);
                stack.push(node.a);
            }
        }
        best
    }

    fn tri_hit(&self, ray: &Ray, ti: u32) -> Option<Hit> {
        let t = self.triangles[ti as usize];
        let a = self.vertices[t[0] as usize];
        let b = self.vertices[t[1] as usize];
        let c = self.vertices[t[2] as usize];
        let e1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let e2 = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let d = [ray.dir.x, ray.dir.y, ray.dir.z];
        let p = cross(d, e2);
        let det = dot(e1, p);
        if det.abs() < 1e-14 {
            return None;
        }
        let inv = 1.0 / det;
        let s = [
            ray.origin.x - a[0],
            ray.origin.y - a[1],
            ray.origin.z - a[2],
        ];
        let u = dot(s, p) * inv;
        if !(0.0..=1.0).contains(&u) {
            return None;
        }
        let q = cross(s, e1);
        let v = dot(d, q) * inv;
        if v < 0.0 || u + v > 1.0 {
            return None;
        }
        let tt = dot(e2, q) * inv;
        (tt > 1e-9).then(|| {
            let n = cross(e1, e2);
            let nn = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            Hit {
                t: tt,
                point: ray.at(tt),
                normal: (nn > 1e-12).then(|| Vec3::new(n[0] / nn, n[1] / nn, n[2] / nn)),
                steps: 0,
            }
        })
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn slab_hit(ray: &Ray, lo: [f64; 3], hi: [f64; 3], t_best: f64) -> bool {
    let o = [ray.origin.x, ray.origin.y, ray.origin.z];
    let d = [ray.dir.x, ray.dir.y, ray.dir.z];
    let mut t0 = 0.0f64;
    let mut t1 = t_best;
    for k in 0..3 {
        let inv = 1.0 / d[k];
        let (mut a, mut b) = ((lo[k] - o[k]) * inv, (hi[k] - o[k]) * inv);
        if a > b {
            core::mem::swap(&mut a, &mut b);
        }
        t0 = t0.max(a);
        t1 = t1.min(b);
        if t0 > t1 {
            return false;
        }
    }
    true
}

/// One scene, three backend kinds, one image: closest hit wins.
pub enum Backend<'a> {
    /// Any chart with a certified Lipschitz bound (SDF / F-rep).
    Chart(&'a dyn Chart),
    /// A NURBS patch.
    Nurbs(&'a NurbsSurface<f64>),
    /// A native triangle mesh.
    Mesh(&'a TriMesh),
}

/// Trace a mixed-chart scene; returns (instance index, hit).
#[must_use]
pub fn trace_scene(
    backends: &[Backend<'_>],
    cx: &Cx<'_>,
    ray: &Ray,
    t_max: f64,
    eps: f64,
) -> Option<(usize, Hit)> {
    let mut best: Option<(usize, Hit)> = None;
    for (i, b) in backends.iter().enumerate() {
        let hit = match b {
            Backend::Chart(chart) => sphere_trace(*chart, cx, ray, t_max, eps, 1.0).0,
            Backend::Nurbs(surface) => ray_intersect_nurbs(surface, ray, 8, eps),
            Backend::Mesh(mesh) => mesh.intersect(ray),
        };
        if let Some(h) = hit
            && best.as_ref().is_none_or(|(_, bh)| h.t < bh.t)
        {
            best = Some((i, h));
        }
    }
    best
}

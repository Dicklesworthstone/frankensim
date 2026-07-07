//! Surface remeshing (plan §7.5): the Botsch–Kobbelt incremental loop —
//! split / collapse / flip / tangential-smooth — measured in a
//! RIEMANNIAN METRIC, with projection back to the source chart.
//! Isotropic remeshing IS the identity-metric special case
//! ([`UniformMetric`]); anisotropic remeshing feeds the same four ops a
//! spatially varying SPD tensor (ultimately FLUX's DWR error metric —
//! stretched elements aligned with solution anisotropy are the
//! difference between 10⁶ and 10⁸ elements at equal accuracy).
//!
//! Architecture: each pass is FUNCTIONAL — connectivity maps
//! (`BTreeMap`/`BTreeSet`) are rebuilt from the triangle list, ops are
//! scheduled in canonical order, and mesh validity is re-established by
//! construction (the conformance battery additionally round-trips
//! through `HalfEdgeMesh::from_triangles` after every round). This
//! trades peak throughput for auditability and P2 determinism — the
//! right trade until the perf-lane bead profiles it.
//!
//! Feature preservation: edges whose dihedral angle exceeds
//! `crease_angle` (plus boundary and non-manifold edges) are LOCKED —
//! never flipped or collapsed; their endpoints never smooth. Split
//! midpoints always chart-project — a no-op on straight creases (the
//! chord lies on the chart), a documented rounding on curved ones.
//! The unit-mesh convention: split above 4/3, collapse below 4/5
//! (metric lengths).

use crate::delaunay::MeshError;
use fs_exec::Cx;
use fs_geom::{Chart, Point3, Vec3};
use fs_rep_mesh::Soup;
use std::collections::{BTreeMap, BTreeSet};

/// A spatially varying SPD metric tensor (row-major, symmetric).
pub trait MetricField: Send + Sync {
    /// The metric at `p`.
    fn metric(&self, p: Point3) -> [[f64; 3]; 3];
}

/// The isotropic case: `M = I / target²` (unit metric length = one
/// target edge length).
#[derive(Debug, Clone, Copy)]
pub struct UniformMetric {
    /// Target edge length (> 0).
    pub target: f64,
}

impl MetricField for UniformMetric {
    fn metric(&self, _p: Point3) -> [[f64; 3]; 3] {
        let s = 1.0 / (self.target * self.target);
        [[s, 0.0, 0.0], [0.0, s, 0.0], [0.0, 0.0, s]]
    }
}

/// Remeshing policy.
#[derive(Debug, Clone, Copy)]
pub struct RemeshOptions {
    /// Full split/collapse/flip/smooth rounds.
    pub iterations: u32,
    /// Dihedral angle (radians) above which an edge is a locked crease.
    pub crease_angle: f64,
    /// Tangential smoothing step (0 disables).
    pub smoothing: f64,
}

impl Default for RemeshOptions {
    fn default() -> Self {
        RemeshOptions {
            iterations: 10,
            crease_angle: 0.7, // ≈ 40°
            smoothing: 0.5,
        }
    }
}

/// Per-run ledger evidence.
#[derive(Debug, Clone, Copy, Default)]
pub struct RemeshStats {
    /// Edge splits applied.
    pub splits: u64,
    /// Edge collapses applied.
    pub collapses: u64,
    /// Edge flips applied.
    pub flips: u64,
    /// Vertex smoothing moves applied.
    pub smooths: u64,
    /// Worst |signed distance| drift from the chart after the run
    /// (0 when no chart is supplied).
    pub worst_chart_drift: f64,
}

impl RemeshStats {
    /// Canonical JSON object.
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"splits\":{},\"collapses\":{},\"flips\":{},\"smooths\":{},\
             \"worst_chart_drift\":{:.3e}}}",
            self.splits, self.collapses, self.flips, self.smooths, self.worst_chart_drift
        )
    }
}

fn sub(a: Point3, b: Point3) -> Vec3 {
    a.delta_from(b)
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}

fn face_normal(p: &[Point3], f: [u32; 3]) -> Vec3 {
    cross(
        sub(p[f[1] as usize], p[f[0] as usize]),
        sub(p[f[2] as usize], p[f[0] as usize]),
    )
}

fn midpoint(a: Point3, b: Point3) -> Point3 {
    Point3::new(
        f64::midpoint(a.x, b.x),
        f64::midpoint(a.y, b.y),
        f64::midpoint(a.z, b.z),
    )
}

/// Metric edge length: `√(eᵀ M(mid) e)`.
fn metric_len(metric: &dyn MetricField, a: Point3, b: Point3) -> f64 {
    let m = metric.metric(midpoint(a, b));
    let e = sub(b, a);
    let v = [e.x, e.y, e.z];
    let mut q = 0.0;
    for (row, &vi) in m.iter().zip(&v) {
        q += vi * (row[0] * v[0] + row[1] * v[1] + row[2] * v[2]);
    }
    q.max(0.0).sqrt()
}

fn ekey(a: u32, b: u32) -> (u32, u32) {
    if a < b { (a, b) } else { (b, a) }
}

/// Newton projection onto the chart's zero set (needs a gradient claim;
/// silently keeps the point where the chart offers none).
fn project(chart: &dyn Chart, mut p: Point3, cx: &Cx<'_>) -> Point3 {
    for _ in 0..3 {
        let s = chart.eval(p, cx);
        let Some(g) = s.gradient else { break };
        p = p.offset(g.scale(-s.signed_distance));
    }
    p
}

struct Passes<'a> {
    positions: Vec<Point3>,
    faces: Vec<[u32; 3]>,
    metric: &'a dyn MetricField,
    chart: Option<&'a dyn Chart>,
    crease_cos: f64,
    stats: RemeshStats,
}

impl Passes<'_> {
    fn edge_faces(&self) -> BTreeMap<(u32, u32), Vec<usize>> {
        let mut map: BTreeMap<(u32, u32), Vec<usize>> = BTreeMap::new();
        for (fi, f) in self.faces.iter().enumerate() {
            for c in 0..3 {
                map.entry(ekey(f[c], f[(c + 1) % 3])).or_default().push(fi);
            }
        }
        map
    }

    /// Locked edges: dihedral creases, boundaries, non-manifold fins.
    fn crease_edges(&self, edge_faces: &BTreeMap<(u32, u32), Vec<usize>>) -> BTreeSet<(u32, u32)> {
        let mut creases = BTreeSet::new();
        for (&e, fs) in edge_faces {
            match fs.as_slice() {
                [f1, f2] => {
                    let n1 = face_normal(&self.positions, self.faces[*f1]);
                    let n2 = face_normal(&self.positions, self.faces[*f2]);
                    let denom = n1.norm() * n2.norm();
                    if denom < 1e-300 || n1.dot(n2) / denom < self.crease_cos {
                        creases.insert(e);
                    }
                }
                _ => {
                    creases.insert(e);
                }
            }
        }
        creases
    }

    fn feature_verts(&self, creases: &BTreeSet<(u32, u32)>) -> BTreeSet<u32> {
        creases.iter().flat_map(|&(a, b)| [a, b]).collect()
    }

    fn len(&self, a: u32, b: u32) -> f64 {
        metric_len(
            self.metric,
            self.positions[a as usize],
            self.positions[b as usize],
        )
    }

    /// Split every edge longer than 4/3 (canonical subdivision patterns,
    /// orientation-preserving; crease midpoints stay on the chord).
    fn split_pass(&mut self, cx: &Cx<'_>) {
        let edge_faces = self.edge_faces();
        let mut mids: BTreeMap<(u32, u32), u32> = BTreeMap::new();
        for &e in edge_faces.keys() {
            if self.len(e.0, e.1) > 4.0 / 3.0 {
                let mut m = midpoint(
                    self.positions[e.0 as usize],
                    self.positions[e.1 as usize],
                );
                // Project ALL midpoints: on straight creases (the zoo's
                // box edges) the chord lies ON the chart, so projection
                // is a no-op there — while transiently mis-classified
                // "creases" on smooth patches would otherwise bake in
                // chord sagitta. Curved-crease rounding is a documented
                // no-claim.
                if let Some(chart) = self.chart {
                    m = project(chart, m, cx);
                }
                self.positions.push(m);
                mids.insert(e, (self.positions.len() - 1) as u32);
            }
        }
        if mids.is_empty() {
            return;
        }
        self.stats.splits += mids.len() as u64;
        let mut out: Vec<[u32; 3]> = Vec::with_capacity(self.faces.len() * 2);
        for &f in &self.faces {
            // Rotate so the marked pattern is canonical.
            let marked: [Option<u32>; 3] =
                core::array::from_fn(|c| mids.get(&ekey(f[c], f[(c + 1) % 3])).copied());
            let count = marked.iter().flatten().count();
            match count {
                0 => out.push(f),
                3 => {
                    let (mab, mbc, mca) = (
                        marked[0].expect("marked"),
                        marked[1].expect("marked"),
                        marked[2].expect("marked"),
                    );
                    out.extend_from_slice(&[
                        [f[0], mab, mca],
                        [f[1], mbc, mab],
                        [f[2], mca, mbc],
                        [mab, mbc, mca],
                    ]);
                }
                _ => {
                    // Rotate so edge (a,b) is marked; for two marks, so
                    // that (c,a) is the UNMARKED edge.
                    let rot = (0..3)
                        .find(|&r| {
                            marked[r].is_some()
                                && (count == 1 || marked[(r + 1) % 3].is_some())
                        })
                        .expect("a marked rotation exists");
                    let (a, b, c) = (f[rot], f[(rot + 1) % 3], f[(rot + 2) % 3]);
                    let mab = marked[rot].expect("marked");
                    if count == 1 {
                        out.extend_from_slice(&[[a, mab, c], [mab, b, c]]);
                    } else {
                        let mbc = marked[(rot + 1) % 3].expect("marked");
                        out.extend_from_slice(&[[mab, b, mbc], [a, mab, mbc], [a, mbc, c]]);
                    }
                }
            }
        }
        self.faces = out;
    }

    /// Collapse edges shorter than 4/5 (shortest first, canonical
    /// tie-break) under the link condition, feature locks, a no-new-long
    /// -edge guard, and a normal-flip guard.
    #[allow(clippy::too_many_lines)] // one guarded transaction per candidate
    fn collapse_pass(&mut self, cx: &Cx<'_>) {
        let edge_faces = self.edge_faces();
        let creases = self.crease_edges(&edge_faces);
        let features = self.feature_verts(&creases);
        let mut neighbors: BTreeMap<u32, BTreeSet<u32>> = BTreeMap::new();
        let mut vert_faces: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
        for (fi, f) in self.faces.iter().enumerate() {
            for c in 0..3 {
                neighbors.entry(f[c]).or_default().insert(f[(c + 1) % 3]);
                neighbors.entry(f[(c + 1) % 3]).or_default().insert(f[c]);
                vert_faces.entry(f[c]).or_default().push(fi);
            }
        }
        let mut candidates: Vec<(u64, (u32, u32))> = edge_faces
            .iter()
            .filter(|(e, fs)| {
                fs.len() == 2
                    && !creases.contains(e)
                    && !features.contains(&e.0)
                    && !features.contains(&e.1)
            })
            .map(|(&e, _)| (self.len(e.0, e.1).to_bits(), e))
            .filter(|&(bits, _)| f64::from_bits(bits) < 0.8)
            .collect();
        candidates.sort_unstable();
        let mut touched: BTreeSet<u32> = BTreeSet::new();
        let mut dead_faces: BTreeSet<usize> = BTreeSet::new();
        for &(_, (u, v)) in &candidates {
            if touched.contains(&u) || touched.contains(&v) {
                continue;
            }
            let (nu, nv) = (&neighbors[&u], &neighbors[&v]);
            if nu.intersection(nv).count() != 2 {
                continue; // link condition: interior manifold edge only
            }
            let mut m = midpoint(self.positions[u as usize], self.positions[v as usize]);
            if let Some(chart) = self.chart {
                m = project(chart, m, cx);
            }
            // Guard: no new over-long edges from the merged vertex.
            let ring: BTreeSet<u32> = nu.union(nv).copied().collect();
            if ring.iter().any(|&w| {
                w != u && w != v && metric_len(self.metric, m, self.positions[w as usize]) > 4.0 / 3.0
            }) {
                continue;
            }
            // Guard: surviving incident faces must not flip normals.
            let mut flips_normal = false;
            for &fi in vert_faces[&u].iter().chain(&vert_faces[&v]) {
                let f = self.faces[fi];
                if f.contains(&u) && f.contains(&v) {
                    continue; // dies with the collapse
                }
                let before = face_normal(&self.positions, f);
                let after_pos: [Point3; 3] = core::array::from_fn(|c| {
                    if f[c] == u || f[c] == v {
                        m
                    } else {
                        self.positions[f[c] as usize]
                    }
                });
                let after = cross(
                    sub(after_pos[1], after_pos[0]),
                    sub(after_pos[2], after_pos[0]),
                );
                if before.dot(after) <= 0.0 {
                    flips_normal = true;
                    break;
                }
            }
            if flips_normal {
                continue;
            }
            // Apply: v merges into u at m.
            self.positions[u as usize] = m;
            for &fi in &vert_faces[&v] {
                let f = &mut self.faces[fi];
                if f.contains(&u) {
                    dead_faces.insert(fi);
                } else {
                    for slot in f.iter_mut() {
                        if *slot == v {
                            *slot = u;
                        }
                    }
                }
            }
            self.stats.collapses += 1;
            touched.insert(u);
            touched.insert(v);
            touched.extend(ring);
        }
        if !dead_faces.is_empty() || self.stats.collapses > 0 {
            let faces = std::mem::take(&mut self.faces);
            self.faces = faces
                .into_iter()
                .enumerate()
                .filter(|(fi, f)| {
                    !dead_faces.contains(fi) && f[0] != f[1] && f[1] != f[2] && f[0] != f[2]
                })
                .map(|(_, f)| f)
                .collect();
        }
    }

    /// Flip interior edges when it strictly reduces the squared valence
    /// deviation from 6 (live valence/edge bookkeeping, canonical order).
    fn flip_pass(&mut self) {
        let edge_faces = self.edge_faces();
        let creases = self.crease_edges(&edge_faces);
        let mut valence: BTreeMap<u32, i64> = BTreeMap::new();
        let mut edge_set: BTreeSet<(u32, u32)> = BTreeSet::new();
        for (&e, _) in &edge_faces {
            *valence.entry(e.0).or_insert(0) += 1;
            *valence.entry(e.1).or_insert(0) += 1;
            edge_set.insert(e);
        }
        let mut pair: BTreeMap<(u32, u32), [usize; 2]> = edge_faces
            .iter()
            .filter(|(_, fs)| fs.len() == 2)
            .map(|(&e, fs)| (e, [fs[0], fs[1]]))
            .collect();
        let candidates: Vec<(u32, u32)> = pair.keys().copied().collect();
        for e in candidates {
            let Some(&[f1, f2]) = pair.get(&e) else {
                continue;
            };
            if creases.contains(&e) {
                continue;
            }
            let (u, v) = e;
            let apex = |fi: usize| -> Option<u32> {
                self.faces[fi].iter().copied().find(|&w| w != u && w != v)
            };
            let (Some(c), Some(d)) = (apex(f1), apex(f2)) else {
                continue;
            };
            if c == d || edge_set.contains(&ekey(c, d)) {
                continue;
            }
            let val = |w: u32| valence.get(&w).copied().unwrap_or(0);
            let dev = |w: u32, delta: i64| {
                let x = val(w) + delta - 6;
                x * x
            };
            let before = dev(u, 0) + dev(v, 0) + dev(c, 0) + dev(d, 0);
            let after = dev(u, -1) + dev(v, -1) + dev(c, 1) + dev(d, 1);
            if after >= before {
                continue;
            }
            // Winding: if f1 traverses u→v, the flipped pair is
            // (u,d,c)/(v,c,d); the mirrored case swaps roles. Guard:
            // both new triangles keep real area and the summed normal
            // direction (no fold-over).
            let dir = {
                let f = self.faces[f1];
                (0..3).any(|k| f[k] == u && f[(k + 1) % 3] == v)
            };
            let (t1, t2): ([u32; 3], [u32; 3]) = if dir {
                ([u, d, c], [v, c, d])
            } else {
                ([u, c, d], [v, d, c])
            };
            let p = &self.positions;
            let add = |x: Vec3, y: Vec3| Vec3::new(x.x + y.x, x.y + y.y, x.z + y.z);
            let sum_old = add(face_normal(p, self.faces[f1]), face_normal(p, self.faces[f2]));
            let (n1, n2) = (face_normal(p, t1), face_normal(p, t2));
            if n1.norm() < 1e-30 || n2.norm() < 1e-30 || sum_old.dot(add(n1, n2)) <= 0.0 {
                continue;
            }
            // Apply + live bookkeeping.
            self.faces[f1] = t1;
            self.faces[f2] = t2;
            edge_set.remove(&e);
            edge_set.insert(ekey(c, d));
            *valence.entry(u).or_insert(0) -= 1;
            *valence.entry(v).or_insert(0) -= 1;
            *valence.entry(c).or_insert(0) += 1;
            *valence.entry(d).or_insert(0) += 1;
            pair.remove(&e);
            pair.insert(ekey(c, d), [f1, f2]);
            for (edge2, fi_old, fi_new) in [
                (ekey(u, d), f2, f1),
                (ekey(v, c), f1, f2),
            ] {
                if let Some(fs) = pair.get_mut(&edge2) {
                    for slot in fs.iter_mut() {
                        if *slot == fi_old {
                            *slot = fi_new;
                        }
                    }
                }
            }
            self.stats.flips += 1;
        }
    }

    /// Jacobi-style tangential smoothing with chart re-projection
    /// (feature vertices stay put).
    fn smooth_pass(&mut self, lambda: f64, cx: &Cx<'_>) {
        if lambda == 0.0 {
            return;
        }
        let edge_faces = self.edge_faces();
        let creases = self.crease_edges(&edge_faces);
        let features = self.feature_verts(&creases);
        let mut neighbors: BTreeMap<u32, BTreeSet<u32>> = BTreeMap::new();
        for &(a, b) in edge_faces.keys() {
            neighbors.entry(a).or_default().insert(b);
            neighbors.entry(b).or_default().insert(a);
        }
        // Area-weighted vertex normals from the CURRENT positions.
        let mut vnormal: BTreeMap<u32, Vec3> = BTreeMap::new();
        for f in &self.faces {
            let n = face_normal(&self.positions, *f);
            for &v in f {
                let e = vnormal.entry(v).or_insert(Vec3::new(0.0, 0.0, 0.0));
                *e = Vec3::new(e.x + n.x, e.y + n.y, e.z + n.z);
            }
        }
        let mut moved: Vec<(u32, Point3)> = Vec::new();
        for (&v, ns) in &neighbors {
            if features.contains(&v) || ns.is_empty() {
                continue;
            }
            let p = self.positions[v as usize];
            let mut g = Vec3::new(0.0, 0.0, 0.0);
            for &w in ns {
                let d = sub(self.positions[w as usize], p);
                g = Vec3::new(g.x + d.x, g.y + d.y, g.z + d.z);
            }
            let g = g.scale(1.0 / ns.len() as f64);
            let n = vnormal.get(&v).copied().unwrap_or(Vec3::new(0.0, 0.0, 0.0));
            let nn = n.norm();
            let d_t = if nn > 1e-300 {
                let n = n.scale(1.0 / nn);
                let along = n.dot(g);
                Vec3::new(g.x - n.x * along, g.y - n.y * along, g.z - n.z * along)
            } else {
                g
            };
            let mut q = p.offset(d_t.scale(lambda));
            if let Some(chart) = self.chart {
                q = project(chart, q, cx);
            }
            moved.push((v, q));
        }
        self.stats.smooths += moved.len() as u64;
        for (v, q) in moved {
            self.positions[v as usize] = q;
        }
    }
}

/// Remesh `soup` toward unit METRIC edge length (isotropic = supply
/// [`UniformMetric`]), optionally projecting onto `chart`. Returns the
/// compacted result and the op ledger.
///
/// # Errors
/// [`MeshError::Cancelled`] between passes.
pub fn remesh(
    soup: &Soup,
    chart: Option<&dyn Chart>,
    metric: &dyn MetricField,
    opts: RemeshOptions,
    cx: &Cx<'_>,
) -> Result<(Soup, RemeshStats), MeshError> {
    let mut passes = Passes {
        positions: soup.positions.clone(),
        faces: soup.triangles.clone(),
        metric,
        chart,
        crease_cos: opts.crease_angle.cos(),
        stats: RemeshStats::default(),
    };
    for _ in 0..opts.iterations {
        cx.checkpoint()?;
        passes.split_pass(cx);
        passes.collapse_pass(cx);
        passes.flip_pass();
        passes.smooth_pass(opts.smoothing, cx);
    }
    // Compact unused vertices (collapse orphans).
    let mut used: BTreeMap<u32, u32> = BTreeMap::new();
    for f in &passes.faces {
        for &v in f {
            let next = used.len() as u32;
            used.entry(v).or_insert(next);
        }
    }
    let mut positions = vec![Point3::new(0.0, 0.0, 0.0); used.len()];
    for (&old, &new) in &used {
        positions[new as usize] = passes.positions[old as usize];
    }
    let triangles: Vec<[u32; 3]> = passes
        .faces
        .iter()
        .map(|f| f.map(|v| used[&v]))
        .collect();
    let mut stats = passes.stats;
    if let Some(chart) = chart {
        for p in &positions {
            let sd = chart.eval(*p, cx).signed_distance.abs();
            stats.worst_chart_drift = stats.worst_chart_drift.max(sd);
        }
    }
    Ok((
        Soup {
            positions,
            triangles,
        },
        stats,
    ))
}

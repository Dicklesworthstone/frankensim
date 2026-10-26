//! Zero-thickness oriented open screens: the unknown is [p], never two
//! coincident closed faces. In H00^(1/2) the free-edge trace vanishes and the
//! physical jump behaves as sqrt(distance). Mesh grading is a declared axis.
//! The weak hypersingular form and singular triangle moments share geometry.rs;
//! all dynamic kernels share helmholtz::WaveNumber, including damped k.
use super::*;
use crate::helmholtz::{
    MAX_DENSE_PANELS, WaveNumber,
    bc::{Admittance, SurfaceBc},
    green_remainder,
    laplace::DampedWaveNumber,
};
use std::collections::{BTreeMap, BTreeSet};

/// A conforming indexed, oriented open two-manifold. Every component must have
/// a free edge. No closure, thickness or orientation repair is synthesized.
pub struct OpenPanels {
    vertices: Vec<Point>,
    triangles: Vec<[usize; 3]>,
    boundary: Vec<bool>,
}
impl OpenPanels {
    /// Reject closed components, duplicate faces, inconsistent orientation,
    /// nonmanifold edges/vertices, unused vertices and degenerate triangles.
    pub fn new(vertices: Vec<Point>, triangles: Vec<[usize; 3]>) -> Result<Self, HelmholtzError> {
        if vertices.is_empty()
            || triangles.is_empty()
            || vertices
                .iter()
                .flatten()
                .any(|v| !v.is_finite() || v.abs() > 1e9)
        {
            return Err(bad("open screen requires finite indexed triangles"));
        }
        let mut edges = BTreeMap::<(usize, usize), Vec<(usize, usize, usize)>>::new();
        let mut faces = BTreeSet::new();
        let mut incident = vec![Vec::new(); vertices.len()];
        for (f, t) in triangles.iter().enumerate() {
            if t.iter().any(|&i| i >= vertices.len()) {
                return Err(bad("open screen vertex index"));
            }
            let mut key = *t;
            key.sort();
            if key[0] == key[1]
                || key[1] == key[2]
                || !faces.insert(key)
                || norm(cross(
                    sub(vertices[t[1]], vertices[t[0]]),
                    sub(vertices[t[2]], vertices[t[0]]),
                )) <= 0.0
            {
                return Err(bad("degenerate or duplicate open screen face"));
            }
            for a in 0..3 {
                let i = t[a];
                let j = t[(a + 1) % 3];
                edges
                    .entry((i.min(j), i.max(j)))
                    .or_default()
                    .push((i, j, f));
                incident[i].push(f);
            }
        }
        let mut boundary = vec![false; vertices.len()];
        let mut adjacent = vec![Vec::new(); triangles.len()];
        for uses in edges.values() {
            match uses.as_slice() {
                &[(i, j, _)] => {
                    boundary[i] = true;
                    boundary[j] = true;
                }
                &[(i, j, f), (u, v, g)] if i == v && j == u => {
                    adjacent[f].push(g);
                    adjacent[g].push(f);
                }
                _ => {
                    return Err(bad(
                        "nonmanifold or inconsistently oriented open screen edge",
                    ));
                }
            }
        }
        // The link of each vertex is one connected fan, with zero or two
        // boundary spokes. Edge incidence alone misses pinched surfaces.
        for (vertex, fan) in incident.iter().enumerate() {
            if fan.is_empty() {
                return Err(bad("unused open screen vertex"));
            }
            let mut link = BTreeMap::<usize, Vec<usize>>::new();
            for &f in fan {
                let other: Vec<_> = triangles[f]
                    .iter()
                    .copied()
                    .filter(|&v| v != vertex)
                    .collect();
                link.entry(other[0]).or_default().push(other[1]);
                link.entry(other[1]).or_default().push(other[0]);
            }
            let ends = link.values().filter(|v| v.len() == 1).count();
            if link.values().any(|v| v.len() > 2) || ends != if boundary[vertex] { 2 } else { 0 } {
                return Err(bad("nonmanifold open screen vertex link"));
            }
            let start = *link.keys().next().ok_or_else(|| bad("empty vertex link"))?;
            let mut seen = BTreeSet::new();
            let mut stack = vec![start];
            while let Some(v) = stack.pop() {
                if seen.insert(v) {
                    stack.extend(&link[&v]);
                }
            }
            if seen.len() != link.len() {
                return Err(bad("disconnected open screen vertex link"));
            }
        }
        let mut seen = vec![false; triangles.len()];
        for start in 0..triangles.len() {
            if !seen[start] {
                let mut stack = vec![start];
                let mut open = false;
                while let Some(f) = stack.pop() {
                    if !seen[f] {
                        seen[f] = true;
                        open |= triangles[f].iter().any(|&v| boundary[v]);
                        stack.extend(&adjacent[f]);
                    }
                }
                if !open {
                    return Err(bad("closed component passed to open screen solver"));
                }
            }
        }
        Ok(Self {
            vertices,
            triangles,
            boundary,
        })
    }
    pub fn vertices(&self) -> &[Point] {
        &self.vertices
    }
    pub fn triangles(&self) -> &[[usize; 3]] {
        &self.triangles
    }
    pub fn is_free_edge_vertex(&self, i: usize) -> bool {
        self.boundary[i]
    }
    pub fn unknowns(&self) -> usize {
        self.boundary.iter().filter(|&&b| !b).count()
    }
}

/// Independent mesh-edge, triangle-quadrature and angular ladder axes. The
/// caller supplies the mesh at edge_level; that declaration never remeshes or
/// silently truncates geometry. Logical storage/work bounds are not RSS caps.
#[derive(Clone, Copy)]
pub struct OpenOptions {
    pub edge_level: u32,
    pub maximum_edge_level: u32,
    /// 4 or 8 Gauss points per Duffy coordinate.
    pub triangle_order: usize,
    /// 0..=2, the existing antipodal angular rule. It bounds the band:
    /// `PreparedScreen::new` refuses k * diameter above ANGULAR_REACH * 2^r.
    pub angular_refinement: u32,
    pub maximum_unknowns: usize,
    pub maximum_cache_bytes: usize,
    pub maximum_work: u64,
    /// Workers for the pair assembly, the outgoing Gram form and the LU
    /// factor. Results are bitwise identical at every count. One worker keeps
    /// the serial loops and their per-pair poll points; more poll per block.
    pub workers: usize,
}

/// The k * diameter that `geometry::directions(r)` resolves, per 2^r. The
/// real-frequency Im N is a Gram form of plane waves exp(ik s.(x-y)) on that
/// rule, and no two points of the screen are farther apart than its vertex
/// diameter, so `PreparedScreen::new` refuses a band above ANGULAR_REACH * 2^r.
/// Up to that chart the rule's worst error on sin(kR)/(kR) is 6.5e-4, 2.8e-3,
/// 2.9e-2 and 2.1e-2 at r = 0, 1, 2 and 3.
///
/// The check at r = 2: a rigid translating free disc over the baffled piston,
/// in dB on axis / 19.8 degrees off axis, on the tests' edge-graded unit
/// discs disc(nr, nt), prepared at kD 35:
///
/// ```text
/// mesh     kD 25            kD 30            kD 35
/// 16x128   +0.159 / +1.102  +0.253 / +0.242  +1.062 / +1.152
/// 20x160   +0.163 / +1.111  +0.267 / +0.260  +1.083 / +1.137
/// 24x160   +0.164 / +1.115  +0.273 / +0.268  +1.092 / +1.130
/// 28x144   +0.164 / +1.115  +0.276 / +0.271  +1.096 / +1.127
/// ```
///
/// The 1.1 dB is converged to 0.04 dB and is there at kD 25, well inside the
/// chart: it is the free disc's edge diffraction, not the rule, and the test
/// admits 1.5 dB. Unguarded, kD 45 solved 24 to 29 dB high at a backward
/// residual near 1e-15, so past the chart nothing else refuses.
pub const ANGULAR_REACH: f64 = 9.0;

/// The largest distance between two vertices. It bounds every pair of
/// points on the triangles.
fn vertex_diameter(vertices: &[Point]) -> f64 {
    let mut widest = 0.0_f64;
    for (i, &a) in vertices.iter().enumerate() {
        for &b in &vertices[i + 1..] {
            widest = widest.max(norm(sub(a, b)));
        }
    }
    widest
}

/// Facet rows per worker in one parallel block; a block's pair values are
/// held until they are scattered in the serial order.
const ROWS_PER_WORKER: usize = 4;

/// `work(row)` for every row of `rows` on `workers` threads (row `start + t`
/// goes to thread `t mod workers`), returned in row order. Scattering the
/// values in that order reproduces the serial loop's bits.
fn rows_on_workers<T: Send>(
    rows: core::ops::Range<usize>,
    workers: usize,
    work: impl Fn(usize) -> T + Sync,
) -> Vec<T> {
    let (start, count) = (rows.start, rows.len());
    let mut out: Vec<Option<T>> = (0..count).map(|_| None).collect();
    std::thread::scope(|scope| {
        let work = &work;
        let handles: Vec<_> = (0..workers.min(count))
            .map(|t| {
                scope.spawn(move || {
                    (t..count)
                        .step_by(workers)
                        .map(|r| (r, work(start + r)))
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        for handle in handles {
            for (r, value) in handle.join().expect("open screen worker panicked") {
                out[r] = Some(value);
            }
        }
    });
    out.into_iter()
        .map(|v| v.expect("every row is computed"))
        .collect()
}
/// Dense open P1 owner. Over-cap open P1 needs a P1 hypersingular iterative
/// route; the closed P0 .109 route is not mislabeled as that implementation.
pub struct PreparedScreen<'a> {
    surface: &'a OpenPanels,
    facets: Vec<Facet>,
    samples: Vec<Vec<Sample>>,
    dof: Vec<Option<usize>>,
    n: usize,
    options: OpenOptions,
    highest_k: f64,
    static_curl: Vec<f64>,
    static_normal: Vec<f64>,
    bytes: usize,
}
impl<'a> PreparedScreen<'a> {
    pub fn new(
        surface: &'a OpenPanels,
        highest_k: f64,
        options: OpenOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, HelmholtzError> {
        poll(&mut cancelled)?;
        if options.edge_level > options.maximum_edge_level {
            return Err(bad("open edge-grading work cap exceeded"));
        }
        let n = surface.unknowns();
        let nf = surface.triangles.len();
        if n == 0 || n > MAX_DENSE_PANELS || nf > MAX_DENSE_PANELS || n > options.maximum_unknowns {
            return Err(HelmholtzError::WorkCap { panels: n.max(nf) });
        }
        if !highest_k.is_finite()
            || highest_k <= 0.0
            || ![4, 8].contains(&options.triangle_order)
            || options.angular_refinement > 2
            || options.maximum_work == 0
            || options.workers == 0
        {
            return Err(bad(
                "invalid open screen frequency or quadrature/work declaration",
            ));
        }
        let mut next = 0;
        let dof = surface
            .boundary
            .iter()
            .map(|&b| {
                if b {
                    None
                } else {
                    let i = next;
                    next += 1;
                    Some(i)
                }
            })
            .collect();
        let mut facets = Vec::with_capacity(nf);
        for &ids in &surface.triangles {
            poll(&mut cancelled)?;
            let triangle = ids.map(|v| surface.vertices[v]);
            let normal = cross(sub(triangle[1], triangle[0]), sub(triangle[2], triangle[0]));
            let area2 = norm(normal);
            facets.push(Facet::new(
                triangle,
                ids,
                area2 / 2.0,
                normal.map(|v| v / area2),
            )?);
        }
        let size = facets.iter().map(|f| f.area).fold(0.0_f64, f64::max).sqrt();
        let ppw = TAU / highest_k / size;
        // Compare the frequency to the declared boundary directly: an exact
        // f6 input must not refuse after a reciprocal round-trip loses an ulp.
        if highest_k > TAU / (MIN_PANELS_PER_WAVELENGTH * size) {
            return Err(HelmholtzError::TooCoarse {
                panels_per_wavelength: ppw,
            });
        }
        // The band the requested angular rule resolves, read from the options
        // the assembly will use, not from the refinement cap.
        let reach = ANGULAR_REACH * f64::from(1_u32 << options.angular_refinement);
        if highest_k * vertex_diameter(&surface.vertices) > reach {
            return Err(bad("open screen band exceeds its angular-resolution chart"));
        }
        let count = n
            .checked_mul(n)
            .ok_or_else(|| bad("open screen storage overflow"))?;
        let bytes = count
            .checked_mul(16)
            .and_then(|v| {
                v.checked_add(
                    nf * (options.triangle_order.pow(2) * core::mem::size_of::<Sample>()
                        + core::mem::size_of::<Facet>()
                        + core::mem::size_of::<Vec<Sample>>())
                        + surface.vertices.len() * core::mem::size_of::<Option<usize>>(),
                )
            })
            .ok_or_else(|| bad("open screen storage overflow"))?;
        if bytes > options.maximum_cache_bytes {
            return Err(HelmholtzError::CacheBudget {
                required_bytes: bytes,
                maximum_bytes: options.maximum_cache_bytes,
            });
        }
        if work(n, nf, options) > u128::from(options.maximum_work) {
            return Err(bad("open hypersingular frequency work cap exceeded"));
        }
        let samples = facets
            .iter()
            .map(|f| samples(f, options.triangle_order))
            .collect();
        let mut out = Self {
            surface,
            facets,
            samples,
            dof,
            n,
            options,
            highest_k,
            static_curl: zeroed(count, 0.0)?,
            static_normal: zeroed(count, 0.0)?,
            bytes,
        };
        out.admit_work()?;
        let both = |out: &Self, i: usize, j: usize| {
            let forward = out.static_pair(i, j);
            let back = if i == j {
                forward
            } else {
                out.static_pair(j, i)
            };
            (forward, back)
        };
        if options.workers <= 1 {
            for i in 0..nf {
                for j in i..nf {
                    poll(&mut cancelled)?;
                    let (forward, back) = both(&out, i, j);
                    out.add_static(i, j, forward, back);
                }
            }
        } else {
            let block = ROWS_PER_WORKER * options.workers;
            for first in (0..nf).step_by(block) {
                poll(&mut cancelled)?;
                let rows = first..nf.min(first + block);
                let pairs = rows_on_workers(rows.clone(), options.workers, |i| {
                    (i..nf).map(|j| both(&out, i, j)).collect::<Vec<_>>()
                });
                for (i, row) in rows.zip(pairs) {
                    for (j, (forward, back)) in (i..nf).zip(row) {
                        out.add_static(i, j, forward, back);
                    }
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(out)
    }
    fn add_static(&mut self, i: usize, j: usize, forward: Pair, back: Pair) {
        let n = self.n;
        let single = 0.5 * (forward.single + back.single);
        for a in 0..3 {
            for b in 0..3 {
                let moment = 0.5 * (forward.shape[a][b] + back.shape[b][a]);
                let curl = -dot(self.facets[i].curl[a], self.facets[j].curl[b]) * single;
                let normal = dot(self.facets[i].normal, self.facets[j].normal) * moment;
                if let (Some(u), Some(v)) = (
                    self.dof[self.facets[i].ids[a]],
                    self.dof[self.facets[j].ids[b]],
                ) {
                    self.static_curl[u * n + v] += curl;
                    self.static_normal[u * n + v] += normal;
                    if i != j {
                        self.static_curl[v * n + u] += curl;
                        self.static_normal[v * n + u] += normal;
                    }
                }
            }
        }
    }
    pub fn unknowns(&self) -> usize {
        self.n
    }
    pub fn cache_bytes(&self) -> usize {
        self.bytes
    }
    pub fn frequency_work(&self) -> u128 {
        work(self.n, self.facets.len(), self.options)
    }
    fn admit_work(&self) -> Result<(), HelmholtzError> {
        if self.frequency_work() > u128::from(self.options.maximum_work) {
            Err(bad("open hypersingular frequency work cap exceeded"))
        } else {
            Ok(())
        }
    }
    fn static_pair(&self, i: usize, j: usize) -> Pair {
        let mut pair = Pair {
            single: 0.0,
            shape: [[0.0; 3]; 3],
        };
        for x in &self.samples[i] {
            let (s, u, _) = static_moments(x.point, &self.facets[j], i == j);
            pair.single += x.weight * s;
            for a in 0..3 {
                for b in 0..3 {
                    pair.shape[a][b] += x.weight * x.shape[a] * u[b];
                }
            }
        }
        pair
    }
    /// The dynamic remainder's single-layer and shape moments of a facet pair.
    fn dynamic_pair<K: WaveNumber>(
        &self,
        wave: K,
        real_frequency: bool,
        i: usize,
        j: usize,
    ) -> (C64, [[C64; 3]; 3]) {
        let mut single = C64::ZERO;
        let mut shape = [[C64::ZERO; 3]; 3];
        for x in &self.samples[i] {
            for y in &self.samples[j] {
                let mut g =
                    green_remainder(wave, norm(sub(x.point, y.point))).scale(x.weight * y.weight);
                if real_frequency {
                    g.im = 0.0;
                }
                single = single + g;
                for a in 0..3 {
                    for b in 0..3 {
                        shape[a][b] = shape[a][b] + g.scale(x.shape[a] * y.shape[b]);
                    }
                }
            }
        }
        (single, shape)
    }
    fn add_dynamic(
        &self,
        h: &mut [C64],
        k2: C64,
        i: usize,
        j: usize,
        single: C64,
        shape: [[C64; 3]; 3],
    ) {
        let n = self.n;
        for a in 0..3 {
            for b in 0..3 {
                if let (Some(u), Some(v)) = (
                    self.dof[self.facets[i].ids[a]],
                    self.dof[self.facets[j].ids[b]],
                ) {
                    let value = k2
                        * shape[a][b].scale(dot(self.facets[i].normal, self.facets[j].normal))
                        - single.scale(dot(self.facets[i].curl[a], self.facets[j].curl[b]));
                    h[u * n + v] = h[u * n + v] + value;
                    if i != j {
                        h[v * n + u] = h[v * n + u] + value;
                    }
                }
            }
        }
    }
    fn coefficients(&self, k: f64, direction: Point) -> Vec<C64> {
        let mut out = vec![C64::ZERO; self.n];
        for (f, samples) in self.facets.iter().zip(&self.samples) {
            let derivative = C64::new(0.0, -k * dot(f.normal, direction));
            for sample in samples {
                let phase = -k * dot(direction, sample.point);
                let value = C64::new(det::cos(phase), det::sin(phase)) * derivative;
                for a in 0..3 {
                    if let Some(u) = self.dof[f.ids[a]] {
                        out[u] = out[u] + value.scale(sample.weight * sample.shape[a]);
                    }
                }
            }
        }
        out
    }
    fn assemble<K: WaveNumber>(
        &self,
        wave: K,
        real_frequency: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HelmholtzError> {
        let k = wave.value();
        let k2 = k * k;
        let n = self.n;
        let mut h = zeroed(n * n, C64::ZERO)?;
        for i in 0..n * n {
            h[i] = C64::from_re(self.static_curl[i]) + k2.scale(self.static_normal[i]);
        }
        let nf = self.facets.len();
        let workers = self.options.workers;
        if workers <= 1 {
            for i in 0..nf {
                for j in i..nf {
                    poll(cancelled)?;
                    let (single, shape) = self.dynamic_pair(wave, real_frequency, i, j);
                    self.add_dynamic(&mut h, k2, i, j, single, shape);
                }
            }
        } else {
            let block = ROWS_PER_WORKER * workers;
            for first in (0..nf).step_by(block) {
                poll(cancelled)?;
                let rows = first..nf.min(first + block);
                let pairs = rows_on_workers(rows.clone(), workers, |i| {
                    (i..nf)
                        .map(|j| self.dynamic_pair(wave, real_frequency, i, j))
                        .collect::<Vec<_>>()
                });
                for (i, row) in rows.zip(pairs) {
                    for (j, (single, shape)) in (i..nf).zip(row) {
                        self.add_dynamic(&mut h, k2, i, j, single, shape);
                    }
                }
            }
        }
        if real_frequency {
            // Im N is the exact plane-wave identity discretized as a Gram
            // form. No clipping of a resistance or of an eigenvalue occurs.
            if workers <= 1 {
                for (direction, weight) in directions(self.options.angular_refinement) {
                    poll(cancelled)?;
                    let c = self.coefficients(k.re, direction);
                    let scale = k.re * weight / (8.0 * PI * PI);
                    for i in 0..n {
                        for j in 0..n {
                            h[i * n + j].im += scale * (c[i].conj() * c[j]).re;
                        }
                    }
                }
            } else {
                // Each entry still takes its directions in the serial order.
                poll(cancelled)?;
                let terms: Vec<(f64, Vec<C64>)> = directions(self.options.angular_refinement)
                    .into_iter()
                    .map(|(direction, weight)| {
                        (
                            k.re * weight / (8.0 * PI * PI),
                            self.coefficients(k.re, direction),
                        )
                    })
                    .collect();
                let mut rows: Vec<Vec<(usize, &mut [C64])>> =
                    (0..workers).map(|_| Vec::new()).collect();
                for (i, row) in h.chunks_mut(n).enumerate() {
                    rows[i % workers].push((i, row));
                }
                std::thread::scope(|scope| {
                    for share in rows {
                        let terms = &terms;
                        scope.spawn(move || {
                            for (i, row) in share {
                                for (scale, c) in terms {
                                    for j in 0..n {
                                        row[j].im += scale * (c[i].conj() * c[j]).re;
                                    }
                                }
                            }
                        });
                    }
                });
                poll(cancelled)?;
            }
        }
        if h.iter().any(|v| !finite(*v)) {
            return Err(bad("nonfinite open hypersingular operator"));
        }
        Ok(h)
    }
    fn weak_image(&self, law: Admittance<'_>, p: &[C64]) -> Result<Vec<C64>, HelmholtzError> {
        let image = law.image(p)?;
        let mut out = vec![C64::ZERO; self.n];
        match law {
            Admittance::Local { beta, .. } => {
                // Integrate beta(x)*p(x)*test(x); M*diag(beta) would lose
                // passivity for a spatially varying nodal law.
                for (f, samples) in self.facets.iter().zip(&self.samples) {
                    for x in samples {
                        let b = (0..3).fold(C64::ZERO, |v, a| v + beta[f.ids[a]].scale(x.shape[a]));
                        let value =
                            (0..3).fold(C64::ZERO, |v, a| v + p[f.ids[a]].scale(x.shape[a]));
                        for a in 0..3 {
                            if let Some(u) = self.dof[f.ids[a]] {
                                out[u] = out[u] + (b * value).scale(x.weight * x.shape[a]);
                            }
                        }
                    }
                }
            }
            _ => {
                for f in &self.facets {
                    for a in 0..3 {
                        if let Some(u) = self.dof[f.ids[a]] {
                            for b in 0..3 {
                                out[u] = out[u]
                                    + image[f.ids[b]].scale(
                                        f.area * if a == b { 1.0 / 6.0 } else { 1.0 / 12.0 },
                                    );
                            }
                        }
                    }
                }
            }
        }
        Ok(out)
    }
    /// Passive sheet law v=v_drive-Y[p]; beta=rho*c*Y. This does not pretend
    /// to model two independently coated faces (which need a second trace).
    pub fn solve(
        &self,
        k: f64,
        medium: Medium,
        bc: SurfaceBc<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<JumpSolution<'_, 'a>, HelmholtzError> {
        if !k.is_finite() || k <= 0.0 || k > self.highest_k {
            return Err(bad("open screen frequency exceeds prepared band"));
        }
        self.solve_wave(k, true, medium, bc, &mut cancelled)
    }
    /// Damped continuation through the SAME WaveNumber seam, k=i*s/c.
    /// No outgoing-power claim is made at a complex frequency.
    pub fn solve_laplace(
        &self,
        s: C64,
        medium: Medium,
        bc: SurfaceBc<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<JumpSolution<'_, 'a>, HelmholtzError> {
        if !finite(s)
            || s.re < 0.0
            || s.abs() == 0.0
            || !medium.sound_speed.is_finite()
            || medium.sound_speed <= 0.0
        {
            return Err(bad("open Laplace solve requires nonzero Re(s)>=0"));
        }
        let k = C64::new(-s.im / medium.sound_speed, s.re / medium.sound_speed);
        if k.abs() > self.highest_k {
            return Err(bad("open Laplace wave exceeds prepared band"));
        }
        self.solve_wave(
            DampedWaveNumber {
                k,
                phase_origin_m: 0.0,
            },
            false,
            medium,
            bc,
            &mut cancelled,
        )
    }
    /// Several real-frequency drives against ONE assembly and LU at `k`: the
    /// results are bitwise those of `drives.len()` separate [`Self::solve`]
    /// calls with the same `admittance`. Each drive passes `solve`'s own
    /// backward-residual and passive-work checks; any failure refuses the whole
    /// call, with no partial results. At most `max_drives` drives, a declared
    /// ceiling (each result keeps only its nodal jump; the factor is n^2).
    pub fn solve_many(
        &self,
        k: f64,
        medium: Medium,
        admittance: Admittance<'_>,
        drives: &[&[C64]],
        max_drives: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<JumpSolution<'_, 'a>>, HelmholtzError> {
        if !k.is_finite() || k <= 0.0 || k > self.highest_k {
            return Err(bad("open screen frequency exceeds prepared band"));
        }
        if drives.is_empty() || drives.len() > max_drives {
            return Err(bad("open multi-drive solve needs 1..=max_drives drives"));
        }
        poll(&mut cancelled)?;
        self.admit_work()?;
        for drive in drives {
            self.validate_drive(
                medium,
                SurfaceBc {
                    drive_velocity: drive,
                    admittance,
                },
            )?;
        }
        let (h, lu) = self.factor_wave(k, true, admittance, &mut cancelled)?;
        let mut out = Vec::with_capacity(drives.len());
        for drive in drives {
            poll(&mut cancelled)?;
            let bc = SurfaceBc {
                drive_velocity: drive,
                admittance,
            };
            out.push(self.finish_drive(k, true, medium, bc, &h, &lu, &mut cancelled)?);
        }
        Ok(out)
    }
    fn validate_drive(&self, medium: Medium, bc: SurfaceBc<'_>) -> Result<(), HelmholtzError> {
        let nv = self.surface.vertices.len();
        bc.admittance.validate(nv)?;
        if bc.drive_velocity.len() != nv
            || bc.drive_velocity.iter().any(|v| !finite(*v))
            || !medium.density.is_finite()
            || medium.density <= 0.0
            || !medium.sound_speed.is_finite()
            || medium.sound_speed <= 0.0
        {
            return Err(bad(
                "open boundary requires finite nodal drive and positive medium",
            ));
        }
        Ok(())
    }
    /// The assembled operator with its admittance columns, and its LU: the
    /// drive-independent part of a solve, shared by `solve_wave` and
    /// `solve_many` so both carry the same bits.
    fn factor_wave<K: WaveNumber>(
        &self,
        wave: K,
        real_frequency: bool,
        admittance: Admittance<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<C64>, fs_la::eigen_complex::LuComplex), HelmholtzError> {
        let nv = self.surface.vertices.len();
        let n = self.n;
        let mut h = self.assemble(wave, real_frequency, cancelled)?;
        let ik = C64::new(0.0, 1.0) * wave.value();
        if !admittance.is_zero() {
            for j in 0..n {
                poll(cancelled)?;
                let mut unit = vec![C64::ZERO; nv];
                for (v, &dof) in self.dof.iter().enumerate() {
                    if dof == Some(j) {
                        unit[v] = C64::ONE;
                    }
                }
                let col = self.weak_image(admittance, &unit)?;
                for i in 0..n {
                    h[i * n + j] = h[i * n + j] + ik * col[i];
                }
            }
        }
        let lu = fs_la::eigen_complex::with_lu_workers(self.options.workers, || {
            lu_complex_with_cancel(&h, n, &mut *cancelled)
        })
        .map_err(lu_error)?;
        Ok((h, lu))
    }
    /// Back-substitute one drive against a shared factor, with the drive's own
    /// backward-residual and passive-work checks.
    #[allow(clippy::too_many_arguments)]
    fn finish_drive<K: WaveNumber>(
        &self,
        wave: K,
        real_frequency: bool,
        medium: Medium,
        bc: SurfaceBc<'_>,
        h: &[C64],
        lu: &fs_la::eigen_complex::LuComplex,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<JumpSolution<'_, 'a>, HelmholtzError> {
        let n = self.n;
        let ik = C64::new(0.0, 1.0) * wave.value();
        let mut rhs = vec![C64::ZERO; n];
        let mut drive_work = vec![C64::ZERO; n];
        for f in &self.facets {
            for a in 0..3 {
                if let Some(u) = self.dof[f.ids[a]] {
                    for b in 0..3 {
                        drive_work[u] = drive_work[u]
                            + bc.drive_velocity[f.ids[b]]
                                .scale(f.area * if a == b { 1.0 / 6.0 } else { 1.0 / 12.0 });
                    }
                }
            }
        }
        for i in 0..n {
            rhs[i] = ik.scale(medium.density * medium.sound_speed) * drive_work[i];
        }
        let mut x = rhs.clone();
        lu.solve_with_cancel(&mut x, &mut *cancelled)
            .map_err(lu_error)?;
        let mut residual = 0.0_f64;
        let mut scale = 0.0_f64;
        for i in 0..n {
            poll(cancelled)?;
            let mut value = C64::ZERO;
            let mut row = rhs[i].abs();
            for j in 0..n {
                value = value + h[i * n + j] * x[j];
                row += h[i * n + j].abs() * x[j].abs();
            }
            residual = residual.max((value - rhs[i]).abs());
            scale = scale.max(row);
        }
        let backward_residual = residual / scale.max(f64::MIN_POSITIVE);
        if !backward_residual.is_finite()
            || backward_residual > 1e-10
            || x.iter().any(|v| !finite(*v))
        {
            return Err(bad("open hypersingular solve backward residual"));
        }
        let jump: Vec<_> = self
            .dof
            .iter()
            .map(|d| d.map_or(C64::ZERO, |i| x[i]))
            .collect();
        let weak = self.weak_image(bc.admittance, &jump)?;
        let absorption = 0.5 / (medium.density * medium.sound_speed)
            * (0..n).map(|i| (x[i].conj() * weak[i]).re).sum::<f64>();
        let driving = 0.5
            * (0..n)
                .map(|i| (x[i] * drive_work[i].conj()).re)
                .sum::<f64>();
        if real_frequency
            && !bc.admittance.declared_active()
            && absorption < -1e-10 * driving.abs().max(f64::MIN_POSITIVE)
        {
            return Err(bad("open nonlocal boundary failed passive-work check"));
        }
        Ok(JumpSolution {
            prepared: self,
            jump,
            k: wave.value(),
            medium,
            real_frequency,
            absorbed_power_w: absorption,
            driving_power_w: driving,
            backward_residual,
            declared_active: bc.admittance.declared_active(),
        })
    }
    fn solve_wave<K: WaveNumber>(
        &self,
        wave: K,
        real_frequency: bool,
        medium: Medium,
        bc: SurfaceBc<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<JumpSolution<'_, 'a>, HelmholtzError> {
        poll(cancelled)?;
        self.admit_work()?;
        self.validate_drive(medium, bc)?;
        let (h, lu) = self.factor_wave(wave, real_frequency, bc.admittance, cancelled)?;
        self.finish_drive(wave, real_frequency, medium, bc, &h, &lu, cancelled)
    }
}

fn work(n: usize, facets: usize, options: OpenOptions) -> u128 {
    let n = n as u128;
    let f = facets as u128;
    let q = options.triangle_order as u128;
    n * n * n
        + f * f * q.pow(4) * 96
        + (64_u128 << (2 * options.angular_refinement)) * (n * n + f * q * q * 24)
}

/// Pressure jump in Pa, with zero free-edge nodal values. A field borrows its
/// prepared geometry, preventing receiver/source cross-wiring.
pub struct JumpSolution<'p, 's> {
    prepared: &'p PreparedScreen<'s>,
    jump: Vec<C64>,
    k: C64,
    medium: Medium,
    real_frequency: bool,
    pub absorbed_power_w: f64,
    pub driving_power_w: f64,
    pub backward_residual: f64,
    pub declared_active: bool,
}
impl JumpSolution<'_, '_> {
    pub fn pressure_jump_pa(&self) -> &[C64] {
        &self.jump
    }
    /// Integral [p] dS [N]; divide by a declared uniform drive [m/s] for Z.
    pub fn force_n(&self) -> C64 {
        self.prepared.facets.iter().fold(C64::ZERO, |sum, f| {
            sum + (self.jump[f.ids[0]] + self.jump[f.ids[1]] + self.jump[f.ids[2]])
                .scale(f.area / 3.0)
        })
    }
    /// r exp(-ikr) p(r*d), in Pa*m. Both sides are represented by D[p].
    pub fn far_field(&self, direction: Point) -> Result<C64, HelmholtzError> {
        if !self.real_frequency
            || direction.iter().any(|v| !v.is_finite())
            || (norm(direction) - 1.0).abs() > 1e-10
        {
            return Err(bad(
                "open far field requires real frequency and unit direction",
            ));
        }
        let c = self.prepared.coefficients(self.k.re, direction);
        Ok(self
            .prepared
            .dof
            .iter()
            .enumerate()
            .fold(C64::ZERO, |sum, (v, d)| {
                sum + d.map_or(C64::ZERO, |u| c[u] * self.jump[v])
            })
            .scale(1.0 / (4.0 * PI)))
    }
    pub fn radiated_power_w(&self) -> Result<f64, HelmholtzError> {
        let mut sum = 0.0;
        for (d, w) in directions(self.prepared.options.angular_refinement) {
            sum += w
                * (self.far_field(d)?.abs().powi(2) + self.far_field(d.map(|v| -v))?.abs().powi(2));
        }
        Ok(sum / (2.0 * self.medium.density * self.medium.sound_speed))
    }
    /// Existing adaptive near-field P1 quadrature; a point on either face or
    /// edge refuses at the declared geometric separation.
    pub fn pressure_at(
        &self,
        points: &[Point],
        separation: f64,
        options: crate::near_field::Options,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HelmholtzError> {
        if !self.real_frequency {
            return Err(bad("open adaptive receivers require real frequency"));
        }
        let triangles: Vec<_> = self.prepared.facets.iter().map(|f| f.triangle).collect();
        let normals: Vec<_> = self.prepared.facets.iter().map(|f| f.normal).collect();
        let jumps: Vec<_> = self
            .prepared
            .facets
            .iter()
            .map(|f| f.ids.map(|v| self.jump[v]))
            .collect();
        crate::near_field::linear::open_pressure(
            self.k.re, &triangles, &normals, &jumps, points, separation, options, cancelled,
        )
    }
}

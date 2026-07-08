//! The Delaunay kernel: BRIO-ordered incremental Bowyer–Watson on exact
//! predicates, with ghost tets carrying the hull.
//!
//! Conventions (load-bearing — the audits pin them):
//! - Every REAL tet `[v0,v1,v2,v3]` is positively oriented:
//!   `orient3d(v0,v1,v2,v3) = Positive` (Shewchuk's convention: v3 below
//!   the CCW-from-above plane v0,v1,v2).
//! - `FACET[i]` lists the facet opposite slot `i`, ordered so the
//!   remaining vertex is Positive ("interior below the facet").
//! - GHOST tets are stored `[g0,g1,g2,GHOST]` where `(g0,g1,g2)` is the
//!   hull facet ordered so a point STRICTLY OUTSIDE the hull across it
//!   is Positive — so the ghost conflict test is one `orient3d_sos`.
//! - Cospherical ties: `insphere = Zero` means NOT in conflict — a
//!   deterministic choice (P2) that keeps the empty-sphere property in
//!   its weak (non-strict) form, which is exactly what the audit checks.
//! - Cavity repair GROWS: a boundary facet not STRICTLY visible from
//!   the new point (plain `orient3d ≠ Positive`) absorbs its outside
//!   neighbor, so no flat tet is ever created and `insphere`'s
//!   positive-orientation precondition always holds.

use fs_exec::{Cancelled, Cx};
use fs_geom::Point3;
use fs_ivl::{Sign, insphere, orient3d, orient3d_sos};
use fs_rep_mesh::{Soup, TetComplex};
use std::collections::BTreeMap;

/// The vertex-at-infinity sentinel carried by hull (ghost) tets.
pub const GHOST: u32 = u32::MAX;

/// Facet opposite slot `i`, ordered so the remaining vertex is Positive.
const FACET: [[usize; 3]; 4] = [[1, 3, 2], [0, 2, 3], [0, 3, 1], [0, 1, 2]];

/// Teaching errors for the meshing entry points.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MeshError {
    /// Fewer than 4 input points.
    TooFewPoints {
        /// How many arrived.
        got: usize,
    },
    /// Every input point lies on one plane (or worse): 3D Delaunay is
    /// undefined — triangulate in 2D instead.
    DegenerateInput,
    /// Cooperative cancellation observed between insertions.
    Cancelled,
}

impl core::fmt::Display for MeshError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            MeshError::TooFewPoints { got } => {
                write!(f, "tetrahedralization needs at least 4 points, got {got}")
            }
            MeshError::DegenerateInput => write!(
                f,
                "all input points are coplanar (exact orient3d test): a 3D \
                 Delaunay tetrahedralization does not exist — use a 2D \
                 triangulation of the common plane instead"
            ),
            MeshError::Cancelled => write!(f, "cancelled between insertions"),
        }
    }
}

impl std::error::Error for MeshError {}

impl From<Cancelled> for MeshError {
    fn from(_: Cancelled) -> Self {
        MeshError::Cancelled
    }
}

/// Build statistics (ledger evidence).
#[derive(Debug, Clone, Copy, Default)]
pub struct DelaunayStats {
    /// Input points.
    pub points_in: u64,
    /// Bitwise-duplicate points skipped (with which receipt).
    pub duplicates_skipped: u64,
    /// Total visibility-walk steps across all insertions.
    pub walk_steps: u64,
    /// Total tets deleted by cavities.
    pub cavity_tets: u64,
    /// Cavity growth repairs (degenerate visibility absorbed).
    pub growth_repairs: u64,
    /// Walks that fell back to exhaustive scan (robustness net).
    pub exhaustive_locates: u64,
    /// Alive real tets at the end.
    pub tets_final: u64,
}

impl DelaunayStats {
    /// Canonical JSON object.
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"points_in\":{},\"duplicates_skipped\":{},\"walk_steps\":{},\
             \"cavity_tets\":{},\"growth_repairs\":{},\"exhaustive_locates\":{},\
             \"tets_final\":{}}}",
            self.points_in,
            self.duplicates_skipped,
            self.walk_steps,
            self.cavity_tets,
            self.growth_repairs,
            self.exhaustive_locates,
            self.tets_final
        )
    }
}

/// What the exact self-audit found (empty = certified clean).
#[derive(Debug, Clone, Default)]
pub struct AuditReport {
    /// Human-readable violations (empty on a healthy mesh).
    pub violations: Vec<String>,
}

impl AuditReport {
    /// True when nothing failed.
    #[must_use]
    pub fn clean(&self) -> bool {
        self.violations.is_empty()
    }
}

#[derive(Debug)]
pub(crate) struct Mesh {
    pub points: Vec<[f64; 3]>,
    /// Tet vertex slots (GHOST allowed in slot 3 only).
    pub tets: Vec<[u32; 4]>,
    /// `adj[t][i]` = tet across the facet opposite slot `i`.
    pub adj: Vec<[u32; 4]>,
    pub alive: Vec<bool>,
    free: Vec<u32>,
    /// Epoch-stamped conflict marks (avoids clearing per insertion).
    mark: Vec<u32>,
    epoch: u32,
    hint: u32,
    pub stats: DelaunayStats,
}

impl Mesh {
    fn coords(&self, v: u32) -> [f64; 3] {
        self.points[v as usize]
    }

    pub(crate) fn is_ghost(&self, t: u32) -> bool {
        self.tets[t as usize][3] == GHOST
    }

    fn alloc(&mut self, verts: [u32; 4]) -> u32 {
        debug_assert!(verts[0] != GHOST && verts[1] != GHOST && verts[2] != GHOST);
        if let Some(t) = self.free.pop() {
            self.tets[t as usize] = verts;
            self.adj[t as usize] = [GHOST; 4];
            self.alive[t as usize] = true;
            self.mark[t as usize] = 0;
            t
        } else {
            self.tets.push(verts);
            self.adj.push([GHOST; 4]);
            self.alive.push(true);
            self.mark.push(0);
            (self.tets.len() - 1) as u32
        }
    }

    fn kill(&mut self, t: u32) {
        self.alive[t as usize] = false;
        self.free.push(t);
    }

    /// SoS orientation of the plane through the facet's three REAL
    /// vertices against query `p` (facet order as given).
    pub(crate) fn facet_sees_sos(&self, f: [u32; 3], p: [f64; 3], p_idx: u32) -> Sign {
        orient3d_sos(
            self.coords(f[0]),
            self.coords(f[1]),
            self.coords(f[2]),
            p,
            [
                u64::from(f[0]),
                u64::from(f[1]),
                u64::from(f[2]),
                u64::from(p_idx),
            ],
        )
    }

    pub(crate) fn facet_verts(&self, t: u32, i: usize) -> [u32; 3] {
        let tv = self.tets[t as usize];
        [tv[FACET[i][0]], tv[FACET[i][1]], tv[FACET[i][2]]]
    }

    /// Is `p` strictly outside the hull facet of ghost `t` — or, when
    /// exactly coplanar, strictly inside the facet's circumcircle IN the
    /// plane (the standard rule: the ghost's "circumsphere" is the
    /// halfspace closed by that disk)? Fully geometric and exact — no
    /// SoS here, so the conflict region is canonical and the cavity
    /// argument (real boundary facets strictly visible) holds.
    fn ghost_conflict(&self, f: [u32; 3], p: [f64; 3]) -> bool {
        let (a, b, c) = (self.coords(f[0]), self.coords(f[1]), self.coords(f[2]));
        match orient3d(a, b, c, p) {
            Sign::Positive => true,
            Sign::Negative => false,
            Sign::Zero => {
                // Project along a fixed axis order; the first
                // non-degenerate projection decides (exact 2D ladder).
                for drop in [2usize, 0, 1] {
                    let q = |v: [f64; 3]| -> [f64; 2] {
                        match drop {
                            2 => [v[0], v[1]],
                            0 => [v[1], v[2]],
                            _ => [v[2], v[0]],
                        }
                    };
                    return match fs_ivl::orient2d(q(a), q(b), q(c)) {
                        Sign::Zero => continue,
                        Sign::Positive => {
                            fs_ivl::incircle(q(a), q(b), q(c), q(p)) == Sign::Positive
                        }
                        Sign::Negative => {
                            fs_ivl::incircle(q(a), q(c), q(b), q(p)) == Sign::Positive
                        }
                    };
                }
                false // degenerate hull facet cannot exist
            }
        }
    }

    /// Is `p` in conflict with tet `t` (strictly inside its
    /// circumsphere; ghosts per [`Mesh::ghost_conflict`])?
    pub(crate) fn in_conflict(&self, t: u32, p: [f64; 3], _p_idx: u32) -> bool {
        let tv = self.tets[t as usize];
        if tv[3] == GHOST {
            self.ghost_conflict([tv[0], tv[1], tv[2]], p)
        } else {
            insphere(
                self.coords(tv[0]),
                self.coords(tv[1]),
                self.coords(tv[2]),
                self.coords(tv[3]),
                p,
            ) == Sign::Positive
        }
    }

    /// Visibility walk from the hint to a tet in conflict with `p`.
    /// Returns the conflict seed. Falls back to an exhaustive scan if
    /// the walk exceeds its budget (logged, deterministic).
    pub(crate) fn locate(&mut self, p: [f64; 3], p_idx: u32) -> u32 {
        let mut t = self.hint;
        if !self.alive[t as usize] || self.is_ghost(t) {
            t = (0..self.tets.len() as u32)
                .find(|&c| self.alive[c as usize] && !self.is_ghost(c))
                .expect("a live real tet always exists");
        }
        let budget = 4 * self.tets.len() as u64 + 64;
        let mut steps = 0u64;
        'walk: loop {
            steps += 1;
            if steps > budget {
                self.stats.exhaustive_locates += 1;
                self.stats.walk_steps += steps;
                return (0..self.tets.len() as u32)
                    .find(|&c| self.alive[c as usize] && self.in_conflict(c, p, p_idx))
                    .expect("a distinct point conflicts with some tet");
            }
            for i in 0..4 {
                let f = self.facet_verts(t, i);
                if self.facet_sees_sos(f, p, p_idx) == Sign::Negative {
                    let n = self.adj[t as usize][i];
                    if self.is_ghost(n) {
                        // Walk exits the hull: the ghost is the seed.
                        self.stats.walk_steps += steps;
                        return n;
                    }
                    t = n;
                    continue 'walk;
                }
            }
            // No facet separates: p is inside t (SoS is total).
            self.stats.walk_steps += steps;
            return t;
        }
    }

    /// Insert point `p_idx` (already appended to `points`). Returns
    /// false when it bitwise-duplicates an existing vertex.
    #[allow(clippy::too_many_lines)] // flood, repair, and rewire are one transaction
    #[allow(clippy::float_cmp)] // duplicate detection is DELIBERATELY bitwise
    pub(crate) fn insert(&mut self, p_idx: u32) -> bool {
        let p = self.coords(p_idx);
        let mut seed = self.locate(p, p_idx);
        // Duplicate guard: identical bits to a vertex of the seed tet
        // (the SoS walk parks a coincident point next to its twin).
        for &v in &self.tets[seed as usize] {
            if v != GHOST && self.coords(v) == p {
                self.stats.duplicates_skipped += 1;
                return false;
            }
        }
        // The SoS walk can exit through a facet p is merely COPLANAR
        // with; the geometric conflict rule may disown that ghost. A
        // strictly-conflicting tet exists for any distinct point (a
        // violated halfspace for outside points, the containing tet
        // otherwise) — find it exhaustively in that rare case.
        if !self.in_conflict(seed, p, p_idx) {
            self.stats.exhaustive_locates += 1;
            seed = (0..self.tets.len() as u32)
                .find(|&c| self.alive[c as usize] && self.in_conflict(c, p, p_idx))
                .expect("a distinct point conflicts with some tet");
        }
        // Conflict flood.
        self.epoch += 1;
        let epoch = self.epoch;
        let mut cavity: Vec<u32> = vec![seed];
        self.mark[seed as usize] = epoch;
        let mut scan = 0;
        while scan < cavity.len() {
            let t = cavity[scan];
            scan += 1;
            for i in 0..4 {
                let n = self.adj[t as usize][i];
                if self.mark[n as usize] != epoch && self.in_conflict(n, p, p_idx) {
                    self.mark[n as usize] = epoch;
                    cavity.push(n);
                }
            }
        }
        // Growth repair: every REAL boundary facet must be STRICTLY
        // visible from p (plain orient3d Positive) so every new tet has
        // positive volume; otherwise absorb the outside neighbor.
        loop {
            let mut grew = false;
            for ci in 0..cavity.len() {
                let t = cavity[ci];
                for i in 0..4 {
                    let n = self.adj[t as usize][i];
                    if self.mark[n as usize] == epoch {
                        continue;
                    }
                    let f = self.facet_verts(t, i);
                    if f.contains(&GHOST) {
                        continue; // ghost side facets are at infinity
                    }
                    let vis = orient3d(self.coords(f[0]), self.coords(f[1]), self.coords(f[2]), p);
                    if vis != Sign::Positive {
                        self.mark[n as usize] = epoch;
                        cavity.push(n);
                        self.stats.growth_repairs += 1;
                        grew = true;
                    }
                }
            }
            if !grew {
                break;
            }
        }
        self.stats.cavity_tets += cavity.len() as u64;
        // Collect boundary facets: (outside neighbor, facet verts as
        // seen from the CAVITY side, i.e. p strictly Positive).
        let mut boundary: Vec<(u32, [u32; 3])> = Vec::new();
        for &t in &cavity {
            for i in 0..4 {
                let n = self.adj[t as usize][i];
                if self.mark[n as usize] != epoch {
                    boundary.push((n, self.facet_verts(t, i)));
                }
            }
        }
        for &t in &cavity {
            self.kill(t);
        }
        // Create one new tet per boundary facet.
        let mut facet_map: BTreeMap<[u32; 3], (u32, usize)> = BTreeMap::new();
        let mut created: Vec<u32> = Vec::with_capacity(boundary.len());
        for &(outside, f) in &boundary {
            let verts = normalize_with(f, p_idx);
            let t = self.alloc(verts);
            created.push(t);
            // Wire the exterior facet: the slot opposite the new vertex.
            let p_slot = self.tets[t as usize]
                .iter()
                .position(|&v| v == p_idx)
                .expect("new vertex present");
            self.adj[t as usize][p_slot] = outside;
            let of = sorted3(f);
            let o_slot = (0..4)
                .find(|&j| sorted3(self.facet_verts(outside, j)) == of)
                .expect("outside neighbor shares the boundary facet");
            self.adj[outside as usize][o_slot] = t;
            // Interior facets (those containing p) pair up via the map.
            for i in 0..4 {
                if i == p_slot {
                    continue;
                }
                let key = sorted3(self.facet_verts(t, i));
                if let Some((u, j)) = facet_map.remove(&key) {
                    self.adj[t as usize][i] = u;
                    self.adj[u as usize][j] = t;
                } else {
                    facet_map.insert(key, (t, i));
                }
            }
        }
        debug_assert!(facet_map.is_empty(), "cavity boundary must close");
        self.hint = *created
            .iter()
            .find(|&&t| !self.is_ghost(t))
            .unwrap_or(&created[0]);
        true
    }
}

/// Sorted facet key (GHOST sorts last).
fn sorted3(f: [u32; 3]) -> [u32; 3] {
    let mut s = f;
    s.sort_unstable();
    s
}

/// Assemble `[f0,f1,f2,p]` and, if GHOST sits among the facet vertices,
/// move it to slot 3 by an EVEN permutation (orientation-preserving).
fn normalize_with(f: [u32; 3], p: u32) -> [u32; 4] {
    let v = [f[0], f[1], f[2], p];
    match v.iter().position(|&x| x == GHOST) {
        None => v,
        Some(0) => [v[1], v[3], v[2], v[0]], // 3-cycle 0→3→1→0
        Some(1) => [v[0], v[2], v[3], v[1]], // 3-cycle 1→3→2→1
        Some(2) => [v[3], v[1], v[0], v[2]], // 3-cycle 2→3→0→2
        Some(_) => unreachable!("new vertex is never GHOST"),
    }
}

/// Morton code (21 bits per axis, interleaved).
fn morton(q: [u32; 3]) -> u64 {
    fn spread(x: u32) -> u64 {
        let mut v = u64::from(x) & 0x1f_ffff;
        v = (v | (v << 32)) & 0x1f_0000_0000_ffff;
        v = (v | (v << 16)) & 0x1f_0000_ff00_00ff;
        v = (v | (v << 8)) & 0x100f_00f0_0f00_f00f;
        v = (v | (v << 4)) & 0x10c3_0c30_c30c_30c3;
        (v | (v << 2)) & 0x1249_2492_4924_9249
    }
    spread(q[0]) | (spread(q[1]) << 1) | (spread(q[2]) << 2)
}

/// BRIO insertion order: deterministic LCG shuffle, then doubling
/// rounds, each round Morton-sorted (spatial locality for the walk
/// without the pathological structure of a single global sweep).
pub(crate) fn brio_order(points: &[[f64; 3]]) -> Vec<u32> {
    let n = points.len();
    let mut order: Vec<u32> = (0..n as u32).collect();
    let mut state = 0x1001_2026_0706_00AAu64;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        state
    };
    for i in (1..n).rev() {
        let j = ((next() >> 32) as usize) % (i + 1);
        order.swap(i, j);
    }
    let (mut lo, mut hi) = (points[0], points[0]);
    for p in points {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let scale: [f64; 3] = core::array::from_fn(|k| {
        let w = hi[k] - lo[k];
        if w > 0.0 {
            f64::from((1u32 << 21) - 1) / w
        } else {
            0.0
        }
    });
    let code = |i: u32| {
        let p = points[i as usize];
        let q: [u32; 3] = core::array::from_fn(|k| {
            ((p[k] - lo[k]) * scale[k]).clamp(0.0, f64::from((1u32 << 21) - 1)) as u32
        });
        morton(q)
    };
    // Rounds: sizes double from a small head to the full set.
    let mut bounds = vec![n];
    let mut m = n;
    while m > 32 {
        m /= 2;
        bounds.push(m);
    }
    bounds.reverse(); // e.g. [32, 64, ..., n]
    let mut start = 0;
    for &b in &bounds {
        order[start..b].sort_by_key(|&i| (code(i), i));
        start = b;
    }
    order
}

/// The finished tetrahedralization: points (input order, Steiner points
/// appended by refinement), live tets, and the exact self-audit.
#[derive(Debug)]
pub struct Tetrahedralization {
    pub(crate) mesh: Mesh,
    /// Index of the first refinement (Steiner) vertex, `== points_in`
    /// until refinement runs.
    pub steiner_from: u32,
}

/// Shared construction prologue: BRIO order, bootstrap quad, first
/// tet + ghosts — used by both the sequential kernel and the colored
/// (parallel) driver so their meshes start bitwise identical.
pub(crate) fn bootstrap_mesh(points: &[Point3]) -> Result<(Mesh, [u32; 4], Vec<u32>), MeshError> {
    if points.len() < 4 {
        return Err(MeshError::TooFewPoints { got: points.len() });
    }
    let pts: Vec<[f64; 3]> = points.iter().map(|p| [p.x, p.y, p.z]).collect();
    let order = brio_order(&pts);
    // Bootstrap: first three affinely independent points in BRIO order,
    // then the first point off their plane (exact tests).
    let mut mesh = Mesh {
        points: pts,
        tets: Vec::new(),
        adj: Vec::new(),
        alive: Vec::new(),
        free: Vec::new(),
        mark: Vec::new(),
        epoch: 0,
        hint: 0,
        stats: DelaunayStats {
            points_in: points.len() as u64,
            ..DelaunayStats::default()
        },
    };
    let quad = bootstrap_quad(&mesh.points, &order).ok_or(MeshError::DegenerateInput)?;
    init_first_tet(&mut mesh, quad);
    Ok((mesh, quad, order))
}

/// Build the Delaunay tetrahedralization of `points` (BRIO order,
/// exact predicates, deterministic). Duplicate points are skipped with
/// a receipt in the stats.
///
/// # Errors
/// [`MeshError::TooFewPoints`], [`MeshError::DegenerateInput`] (exact
/// all-coplanar detection), [`MeshError::Cancelled`].
pub fn delaunay(points: &[Point3], cx: &Cx<'_>) -> Result<Tetrahedralization, MeshError> {
    let (mut mesh, quad, order) = bootstrap_mesh(points)?;
    let mut inserted = 0u64;
    for &i in &order {
        if quad.contains(&i) {
            continue;
        }
        mesh.insert(i);
        inserted += 1;
        if inserted.is_multiple_of(256) {
            cx.checkpoint()?;
        }
    }
    mesh.stats.tets_final = (0..mesh.tets.len() as u32)
        .filter(|&t| mesh.alive[t as usize] && !mesh.is_ghost(t))
        .count() as u64;
    Ok(Tetrahedralization {
        mesh,
        steiner_from: points.len() as u32,
    })
}

/// First 4 points (in BRIO order) that span 3D, by exact tests.
#[allow(clippy::float_cmp)] // distinctness is DELIBERATELY bitwise
pub(crate) fn bootstrap_quad(pts: &[[f64; 3]], order: &[u32]) -> Option<[u32; 4]> {
    let a = order[0];
    // First point not equal to a.
    let b = *order
        .iter()
        .find(|&&i| pts[i as usize] != pts[a as usize])?;
    // First point not collinear with (a, b): some axis-projection pair
    // of orient2d... use orient3d against a probe: simpler exact route —
    // scan for c making a nonzero triangle normal via orient3d with a
    // fourth scan point is circular; instead test collinearity exactly:
    // cross-product expansion sign via orient2d on the three axis
    // projections would need fs-ivl orient2d; cheaper: c is non-collinear
    // iff SOME d has orient3d(a,b,c,d) != Zero. Fuse the scans: find
    // (c, d) with orient3d(a,b,c,d) != Zero.
    let mut c_candidates = order.iter().filter(|&&i| i != a && i != b);
    for &c in c_candidates.by_ref() {
        for &d in order {
            if d == a || d == b || d == c {
                continue;
            }
            if orient3d(
                pts[a as usize],
                pts[b as usize],
                pts[c as usize],
                pts[d as usize],
            ) != Sign::Zero
            {
                return Some([a, b, c, d]);
            }
        }
        // c collinear with every d ⇒ try the next c... but a full inner
        // scan already proves THIS c yields nothing; continue.
    }
    None
}

/// Create the first real tet (swapped Positive) and its 4 ghosts.
pub(crate) fn init_first_tet(mesh: &mut Mesh, quad: [u32; 4]) {
    let [a, b, c, d] = quad;
    let verts = if orient3d(
        mesh.points[a as usize],
        mesh.points[b as usize],
        mesh.points[c as usize],
        mesh.points[d as usize],
    ) == Sign::Positive
    {
        [a, b, c, d]
    } else {
        [b, a, c, d]
    };
    let t = mesh.alloc(verts);
    // Ghost per facet: the facet REVERSED (outside strictly Positive).
    let ghosts: [u32; 4] = core::array::from_fn(|i| {
        let f = mesh.facet_verts(t, i);
        let g = mesh.alloc_ghost([f[0], f[2], f[1]]);
        mesh.adj[t as usize][i] = g;
        mesh.adj[g as usize][3] = t;
        g
    });
    // Wire ghost side facets to each other via the shared-edge map.
    let mut map: BTreeMap<[u32; 3], (u32, usize)> = BTreeMap::new();
    for &g in &ghosts {
        for i in 0..3 {
            let key = sorted3(mesh.facet_verts(g, i));
            if let Some((u, j)) = map.remove(&key) {
                mesh.adj[g as usize][i] = u;
                mesh.adj[u as usize][j] = g;
            } else {
                map.insert(key, (g, i));
            }
        }
    }
    debug_assert!(map.is_empty(), "initial hull must close");
    mesh.hint = t;
}

impl Mesh {
    fn alloc_ghost(&mut self, hull_facet: [u32; 3]) -> u32 {
        if let Some(t) = self.free.pop() {
            self.tets[t as usize] = [hull_facet[0], hull_facet[1], hull_facet[2], GHOST];
            self.adj[t as usize] = [GHOST; 4];
            self.alive[t as usize] = true;
            self.mark[t as usize] = 0;
            t
        } else {
            self.tets
                .push([hull_facet[0], hull_facet[1], hull_facet[2], GHOST]);
            self.adj.push([GHOST; 4]);
            self.alive.push(true);
            self.mark.push(0);
            (self.tets.len() - 1) as u32
        }
    }
}

impl Tetrahedralization {
    /// The vertex coordinates (input order; Steiner points appended).
    #[must_use]
    pub fn points(&self) -> Vec<Point3> {
        self.mesh
            .points
            .iter()
            .map(|p| Point3::new(p[0], p[1], p[2]))
            .collect()
    }

    /// Live real tets, deterministically ordered (sorted by their
    /// sorted vertex tuples), each positively oriented as stored.
    #[must_use]
    pub fn tets(&self) -> Vec<[u32; 4]> {
        let mut out: Vec<[u32; 4]> = (0..self.mesh.tets.len())
            .filter(|&t| self.mesh.alive[t] && self.mesh.tets[t][3] != GHOST)
            .map(|t| self.mesh.tets[t])
            .collect();
        out.sort_unstable_by_key(|t| {
            let mut s = *t;
            s.sort_unstable();
            s
        });
        out
    }

    /// Build statistics.
    #[must_use]
    pub fn stats(&self) -> DelaunayStats {
        self.mesh.stats
    }

    /// The boundary hull as an outward-oriented triangle soup.
    #[must_use]
    pub fn hull(&self) -> Soup {
        let mut triangles = Vec::new();
        for t in 0..self.mesh.tets.len() {
            if self.mesh.alive[t] && self.mesh.tets[t][3] == GHOST {
                let [a, b, c, _] = self.mesh.tets[t];
                // Stored (a,b,c) has outside-Positive = outside "below"
                // (Shewchuk); outward CCW-from-outside is the reverse.
                triangles.push([a, c, b]);
            }
        }
        triangles.sort_unstable();
        Soup {
            positions: self.points(),
            triangles,
        }
    }

    /// The oriented tet complex (δδ = 0 integration point).
    #[must_use]
    pub fn complex(&self) -> TetComplex {
        TetComplex::from_tets(self.mesh.points.len(), self.tets())
    }

    /// Exact self-audit: orientation, mutual adjacency, LOCAL Delaunay
    /// on every internal facet (the Delaunay lemma lifts it to global),
    /// Euler characteristic of the ball, hull closure, and exact hull
    /// convexity. With `full_insphere`, additionally the O(n·t) global
    /// empty-circumsphere check (fixture-scale belt and braces).
    #[must_use]
    #[allow(clippy::too_many_lines)] // one auditing law per block
    pub fn audit(&self, full_insphere: bool) -> AuditReport {
        let m = &self.mesh;
        let mut violations = Vec::new();
        let live: Vec<u32> = (0..m.tets.len() as u32)
            .filter(|&t| m.alive[t as usize])
            .collect();
        // Orientation + adjacency reciprocity + facet agreement.
        for &t in &live {
            let tv = m.tets[t as usize];
            if tv[3] != GHOST {
                let s = orient3d(
                    m.points[tv[0] as usize],
                    m.points[tv[1] as usize],
                    m.points[tv[2] as usize],
                    m.points[tv[3] as usize],
                );
                if s != Sign::Positive {
                    violations.push(format!("tet {t} not positively oriented ({s:?})"));
                }
            }
            for i in 0..4 {
                let n = m.adj[t as usize][i];
                if n == GHOST || !m.alive[n as usize] {
                    violations.push(format!("tet {t} slot {i} adjacency dead/unset"));
                    continue;
                }
                let f = sorted3(m.facet_verts(t, i));
                if !(0..4).any(|j| m.adj[n as usize][j] == t && sorted3(m.facet_verts(n, j)) == f) {
                    violations.push(format!("tets {t}/{n} adjacency not mutual on facet"));
                }
            }
        }
        // Local Delaunay: neighbor apex never STRICTLY inside.
        for &t in &live {
            let tv = m.tets[t as usize];
            if tv[3] == GHOST {
                continue;
            }
            for i in 0..4 {
                let n = m.adj[t as usize][i];
                let nv = m.tets[n as usize];
                if nv[3] == GHOST {
                    continue;
                }
                let f = sorted3(m.facet_verts(t, i));
                let apex = nv
                    .iter()
                    .copied()
                    .find(|v| !f.contains(v))
                    .expect("neighbor has an apex");
                if insphere(
                    m.points[tv[0] as usize],
                    m.points[tv[1] as usize],
                    m.points[tv[2] as usize],
                    m.points[tv[3] as usize],
                    m.points[apex as usize],
                ) == Sign::Positive
                {
                    violations.push(format!(
                        "local Delaunay violated: apex {apex} strictly inside tet {t}"
                    ));
                }
            }
        }
        // Optional global audit.
        if full_insphere {
            for &t in &live {
                let tv = m.tets[t as usize];
                if tv[3] == GHOST {
                    continue;
                }
                for q in 0..m.points.len() as u32 {
                    if tv.contains(&q) {
                        continue;
                    }
                    if insphere(
                        m.points[tv[0] as usize],
                        m.points[tv[1] as usize],
                        m.points[tv[2] as usize],
                        m.points[tv[3] as usize],
                        m.points[q as usize],
                    ) == Sign::Positive
                    {
                        violations.push(format!(
                            "global empty-sphere violated: point {q} inside tet {t}"
                        ));
                    }
                }
            }
        }
        // Euler characteristic of the tetrahedralized ball: V−E+F−T = 1.
        {
            let mut verts = std::collections::BTreeSet::new();
            let mut edges = std::collections::BTreeSet::new();
            let mut faces = std::collections::BTreeSet::new();
            let mut ntet = 0i64;
            for &t in &live {
                let tv = m.tets[t as usize];
                if tv[3] == GHOST {
                    continue;
                }
                ntet += 1;
                for &v in &tv {
                    verts.insert(v);
                }
                for a in 0..4 {
                    for b in (a + 1)..4 {
                        edges.insert((tv[a].min(tv[b]), tv[a].max(tv[b])));
                    }
                }
                for i in 0..4 {
                    faces.insert(sorted3(m.facet_verts(t, i)));
                }
            }
            let count = |n: usize| i64::try_from(n).expect("mesh far below i64::MAX");
            let chi = count(verts.len()) - count(edges.len()) + count(faces.len()) - ntet;
            if chi != 1 {
                violations.push(format!("Euler characteristic {chi} != 1"));
            }
        }
        // Hull: closed (each directed edge once) and exactly convex.
        {
            let mut directed = std::collections::BTreeSet::new();
            for &t in &live {
                let tv = m.tets[t as usize];
                if tv[3] != GHOST {
                    continue;
                }
                let f = [tv[0], tv[1], tv[2]];
                for k in 0..3 {
                    if !directed.insert((f[k], f[(k + 1) % 3])) {
                        violations.push(format!(
                            "hull edge ({},{}) traversed twice",
                            f[k],
                            f[(k + 1) % 3]
                        ));
                    }
                }
                for q in 0..m.points.len() as u32 {
                    if f.contains(&q) {
                        continue;
                    }
                    if orient3d(
                        m.points[f[0] as usize],
                        m.points[f[1] as usize],
                        m.points[f[2] as usize],
                        m.points[q as usize],
                    ) == Sign::Positive
                    {
                        violations.push(format!(
                            "hull not convex: point {q} strictly outside facet of ghost {t}"
                        ));
                    }
                }
            }
            let mut ok = true;
            for &(a, b) in &directed {
                ok &= directed.contains(&(b, a));
            }
            if !ok {
                violations.push("hull surface not closed".to_string());
            }
        }
        AuditReport { violations }
    }
}

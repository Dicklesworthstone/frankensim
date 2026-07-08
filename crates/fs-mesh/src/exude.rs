//! Sliver removal by deterministic Steiner perturbation (bead uee3):
//! slivers — tets whose radius-edge ratio is acceptable but whose
//! minimum DIHEDRAL angle collapses — get their STEINER vertices
//! nudged by small deterministic offsets and the mesh rebuilt through
//! the exact kernel; a round is kept only if the sliver count strictly
//! drops and the exact audit stays clean. Input vertices are NEVER
//! touched (they are data), which is why input-only slivers are
//! counted as protected rather than chased.
//!
//! This is the perturbation flavor of Edelsbrunner-style exudation:
//! the weighted-Delaunay pump assigns weights, this assigns positions
//! — both escape the sliver's measure-zero degeneracy; weights need a
//! weighted exact predicate (recorded no-claim), positions reuse the
//! unweighted kernel and its audit unchanged.

use crate::delaunay::{GHOST, MeshError, Tetrahedralization, delaunay};
use fs_exec::Cx;
use fs_geom::Point3;

/// Exudation policy.
#[derive(Debug, Clone, Copy)]
pub struct ExudeOptions {
    /// Minimum acceptable dihedral angle (degrees).
    pub dihedral_min_deg: f64,
    /// Perturbation rounds.
    pub rounds: u32,
    /// Jitter amplitude as a fraction of the sliver's shortest edge.
    pub jitter: f64,
}

impl Default for ExudeOptions {
    fn default() -> Self {
        ExudeOptions {
            dihedral_min_deg: 5.0,
            rounds: 8,
            jitter: 0.05,
        }
    }
}

/// Exudation evidence.
#[derive(Debug, Clone, Copy, Default)]
pub struct ExudeStats {
    /// Slivers before.
    pub slivers_before: u32,
    /// Slivers after.
    pub slivers_after: u32,
    /// Slivers whose vertices are all INPUT points (protected).
    pub input_protected: u32,
    /// Rounds accepted.
    pub rounds_accepted: u32,
    /// Full rebuilds performed.
    pub rebuilds: u32,
    /// Worst (smallest) dihedral before, degrees.
    pub worst_dihedral_before: f64,
    /// Worst dihedral after, degrees.
    pub worst_dihedral_after: f64,
}

impl ExudeStats {
    /// Canonical JSON row.
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"slivers_before\":{},\"slivers_after\":{},\"input_protected\":{},\
             \"rounds_accepted\":{},\"rebuilds\":{},\"worst_dihedral_before\":{:.3},\
             \"worst_dihedral_after\":{:.3}}}",
            self.slivers_before,
            self.slivers_after,
            self.input_protected,
            self.rounds_accepted,
            self.rebuilds,
            self.worst_dihedral_before,
            self.worst_dihedral_after
        )
    }
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
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

/// Minimum dihedral angle of a tet (degrees).
fn min_dihedral_deg(p: &[[f64; 3]; 4]) -> f64 {
    // Face normals opposite each vertex.
    let n = [
        cross(sub(p[2], p[1]), sub(p[3], p[1])), // opposite 0
        cross(sub(p[3], p[0]), sub(p[2], p[0])), // opposite 1
        cross(sub(p[1], p[0]), sub(p[3], p[0])), // opposite 2
        cross(sub(p[2], p[0]), sub(p[1], p[0])), // opposite 3
    ];
    let mut worst = 180.0f64;
    for i in 0..4 {
        for j in (i + 1)..4 {
            let (a, b) = (n[i], n[j]);
            let c = dot(a, b) / (dot(a, a).sqrt() * dot(b, b).sqrt()).max(1e-300);
            // Dihedral along the shared edge is π − angle(normals).
            let dihedral = 180.0 - c.clamp(-1.0, 1.0).acos().to_degrees();
            worst = worst.min(dihedral);
        }
    }
    worst
}

/// Sliver census: (count, worst dihedral, offending Steiner vertices,
/// input-only slivers).
fn census(tetra: &Tetrahedralization, threshold_deg: f64) -> (u32, f64, Vec<u32>, u32) {
    let pts = tetra.points();
    let mut count = 0u32;
    let mut worst = 180.0f64;
    let mut steiner: Vec<u32> = Vec::new();
    let mut input_only = 0u32;
    for tet in tetra.tets() {
        if tet[3] == GHOST {
            continue;
        }
        let q: [[f64; 3]; 4] = core::array::from_fn(|k| {
            [
                pts[tet[k] as usize].x,
                pts[tet[k] as usize].y,
                pts[tet[k] as usize].z,
            ]
        });
        let d = min_dihedral_deg(&q);
        worst = worst.min(d);
        if d < threshold_deg {
            count += 1;
            let mut any = false;
            for &v in &tet {
                if v >= tetra.steiner_from {
                    steiner.push(v);
                    any = true;
                }
            }
            if !any {
                input_only += 1;
            }
        }
    }
    steiner.sort_unstable();
    steiner.dedup();
    (count, worst, steiner, input_only)
}

/// Deterministic unit jitter direction from (vertex, round).
fn jitter_dir(v: u32, round: u32) -> [f64; 3] {
    let mut h = 0xcbf2_9ce4_8422_2325u64 ^ (u64::from(v) << 17) ^ u64::from(round);
    let mut next = || {
        h ^= h << 13;
        h ^= h >> 7;
        h ^= h << 17;
        #[allow(clippy::cast_precision_loss)]
        {
            (h >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        }
    };
    let d = [next(), next(), next()];
    let n = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt().max(1e-30);
    [d[0] / n, d[1] / n, d[2] / n]
}

/// Run exudation; the tetrahedralization is rebuilt through the exact
/// kernel each accepted round. Input points are untouched.
///
/// # Errors
/// Kernel errors propagate ([`MeshError::Cancelled`] between rounds).
pub fn exude(
    tetra: &mut Tetrahedralization,
    opts: ExudeOptions,
    cx: &Cx<'_>,
) -> Result<ExudeStats, MeshError> {
    let mut stats = ExudeStats::default();
    let (before, worst_before, _, _) = census(tetra, opts.dihedral_min_deg);
    stats.slivers_before = before;
    stats.worst_dihedral_before = worst_before;
    let steiner_from = tetra.steiner_from;
    let mut best_count = before;
    for round in 0..opts.rounds {
        cx.checkpoint()?;
        let (count, _, movers, input_only) = census(tetra, opts.dihedral_min_deg);
        stats.input_protected = input_only;
        if count == 0 || movers.is_empty() {
            break;
        }
        // Perturb the offending Steiner vertices.
        let pts = tetra.points();
        let mut moved: Vec<Point3> = pts.clone();
        for &v in &movers {
            // Local scale: shortest edge among this vertex's uses in
            // sliver tets (conservative: global shortest incident).
            let mut scale = f64::INFINITY;
            for tet in tetra.tets() {
                if tet[3] == GHOST || !tet.contains(&v) {
                    continue;
                }
                for i in 0..4 {
                    for j in (i + 1)..4 {
                        let a = pts[tet[i] as usize];
                        let b = pts[tet[j] as usize];
                        let d = sub([a.x, a.y, a.z], [b.x, b.y, b.z]);
                        scale = scale.min(dot(d, d).sqrt());
                    }
                }
            }
            if !scale.is_finite() {
                continue;
            }
            let dir = jitter_dir(v, round);
            let p = moved[v as usize];
            moved[v as usize] = Point3::new(
                p.x + opts.jitter * scale * dir[0],
                p.y + opts.jitter * scale * dir[1],
                p.z + opts.jitter * scale * dir[2],
            );
        }
        // Rebuild through the exact kernel and re-audit.
        let rebuilt = delaunay(&moved, cx)?;
        stats.rebuilds += 1;
        let mut candidate = rebuilt;
        candidate.steiner_from = steiner_from;
        let audit = candidate.audit(true);
        let (new_count, _, _, _) = census(&candidate, opts.dihedral_min_deg);
        if audit.clean() && new_count < best_count {
            best_count = new_count;
            *tetra = candidate;
            stats.rounds_accepted += 1;
            if new_count == 0 {
                break;
            }
        }
    }
    let (after, worst_after, _, input_only) = census(tetra, opts.dihedral_min_deg);
    stats.slivers_after = after;
    stats.worst_dihedral_after = worst_after;
    stats.input_protected = input_only;
    Ok(stats)
}

//! Deterministic parallel domain coloring (bead uee3 item 4):
//! READ-PARALLEL, APPLY-CANONICAL prefix batches. Each window of BRIO
//! order gets every point's conflict region (cavity + growth repair +
//! one-ring, mirroring the kernel's insert transaction exactly)
//! computed READ-ONLY across scoped threads; points are then colored
//! FLIP-SAFELY (k = 1 + the largest overlapping color): same-color
//! members are pairwise disjoint AND any pair whose application order
//! flips relative to BRIO is disjoint — so cospherical TIE groups
//! (whose weak-Delaunay resolution is order-dependent) keep their
//! original order and the finished mesh merges canonically with the
//! sequential kernel on EVERY input. (Two rejected designs, both
//! measured: first-fit coloring flipped tied pairs and diverged on
//! the 6×6×6 grid; stop-at-first-clash prefix batching preserved raw
//! bitwise order but BRIO locality collapsed its batch width to ~3.)
//! Thread count can only change the wall clock, never a bit of the
//! result; within-color commutativity is gated adversarially by the
//! battery (reversed application, canonical equality).

use crate::delaunay::{GHOST, Mesh, MeshError, Tetrahedralization, bootstrap_mesh};
use fs_exec::Cx;
use std::collections::BTreeSet;

/// Batching ledger.
#[derive(Debug, Clone, Copy, Default)]
pub struct ColoredStats {
    /// Prefix batches applied.
    pub batches: u64,
    /// Largest batch (parallel width evidence).
    pub largest_batch: u64,
    /// Batches of size one (the serial tail of degenerate inputs).
    pub singleton_batches: u64,
    /// Points scheduled (excludes the bootstrap quad).
    pub points: u64,
    /// Thread count used for the read phase.
    pub threads: u64,
}

impl ColoredStats {
    /// Canonical JSON ledger row.
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"batches\":{},\"largest_batch\":{},\"singleton_batches\":{},\
             \"points\":{},\"threads\":{}}}",
            self.batches, self.largest_batch, self.singleton_batches, self.points, self.threads
        )
    }
}

/// Read-only visibility walk to a conflict seed (mirrors
/// `Mesh::locate` without stats/hint mutation).
fn locate_ro(mesh: &Mesh, p: [f64; 3], p_idx: u32) -> u32 {
    let mut t = (0..mesh.tets.len() as u32)
        .find(|&c| mesh.alive[c as usize] && !mesh.is_ghost(c))
        .expect("a live real tet always exists");
    let budget = 4 * mesh.tets.len() as u64 + 64;
    let mut steps = 0u64;
    'walk: loop {
        steps += 1;
        if steps > budget {
            return (0..mesh.tets.len() as u32)
                .find(|&c| mesh.alive[c as usize] && mesh.in_conflict(c, p, p_idx))
                .expect("a distinct point conflicts with some tet");
        }
        for i in 0..4 {
            let f = mesh.facet_verts(t, i);
            if mesh.facet_sees_sos(f, p, p_idx) == fs_ivl::Sign::Negative {
                let n = mesh.adj[t as usize][i];
                if mesh.is_ghost(n) {
                    return n;
                }
                t = n;
                continue 'walk;
            }
        }
        return t;
    }
}

/// The conflict REGION of a point against the current mesh: the
/// cavity (flood + growth repair, mirroring `Mesh::insert`) plus its
/// one-ring of outside neighbors (insertion rewires their adjacency
/// rows, so disjoint REGIONS — not just cavities — is what makes
/// same-color insertions commute). `None` marks a bitwise duplicate
/// of an existing vertex.
#[allow(clippy::float_cmp)] // duplicate detection is DELIBERATELY bitwise
fn conflict_region(mesh: &Mesh, p_idx: u32) -> Option<BTreeSet<u32>> {
    let p = mesh.points[p_idx as usize];
    let mut seed = locate_ro(mesh, p, p_idx);
    for &v in &mesh.tets[seed as usize] {
        if v != GHOST && mesh.points[v as usize] == p {
            return None;
        }
    }
    if !mesh.in_conflict(seed, p, p_idx) {
        seed = (0..mesh.tets.len() as u32)
            .find(|&c| mesh.alive[c as usize] && mesh.in_conflict(c, p, p_idx))
            .expect("a distinct point conflicts with some tet");
    }
    // Flood.
    let mut cavity: Vec<u32> = vec![seed];
    let mut in_cavity: BTreeSet<u32> = BTreeSet::from([seed]);
    let mut scan = 0;
    while scan < cavity.len() {
        let t = cavity[scan];
        scan += 1;
        for i in 0..4 {
            let n = mesh.adj[t as usize][i];
            if !in_cavity.contains(&n) && mesh.in_conflict(n, p, p_idx) {
                in_cavity.insert(n);
                cavity.push(n);
            }
        }
    }
    // Growth repair (identical rule to the kernel's).
    loop {
        let mut grew = false;
        let mut ci = 0;
        while ci < cavity.len() {
            let t = cavity[ci];
            ci += 1;
            for i in 0..4 {
                let n = mesh.adj[t as usize][i];
                if in_cavity.contains(&n) {
                    continue;
                }
                let f = mesh.facet_verts(t, i);
                if f.contains(&GHOST) {
                    continue;
                }
                let vis = fs_ivl::orient3d(
                    mesh.points[f[0] as usize],
                    mesh.points[f[1] as usize],
                    mesh.points[f[2] as usize],
                    p,
                );
                if vis != fs_ivl::Sign::Positive {
                    in_cavity.insert(n);
                    cavity.push(n);
                    grew = true;
                }
            }
        }
        if !grew {
            break;
        }
    }
    // One-ring: outside neighbors (their adjacency rows get rewired).
    let mut region = in_cavity.clone();
    for &t in &cavity {
        for i in 0..4 {
            let n = mesh.adj[t as usize][i];
            if !in_cavity.contains(&n) {
                region.insert(n);
            }
        }
    }
    Some(region)
}

/// Flip-safe color assignment: point p joins color
/// k = 1 + max{ j : region(p) ∩ occupied[j] ≠ ∅ } (0 when disjoint
/// from everything). Two properties BY CONSTRUCTION: (i) same-color
/// members are pairwise region-disjoint (the parallel-apply batch
/// structure, gated by reversed application); (ii) every cross-color
/// pair whose application order flips relative to BRIO is
/// region-disjoint — overlapping points (including cospherical TIE
/// groups, whose regions coincide) cascade through strictly
/// increasing colors and keep their original order, which is why the
/// degenerate-grid gate passes canonically (a first-fit scheduler
/// flipped tied pairs and measurably diverged there). Duplicates
/// (None regions) have an empty footprint and land in color 0.
fn assign_colors(regions: &[(u32, Option<BTreeSet<u32>>)]) -> Vec<Vec<u32>> {
    let mut colors: Vec<(Vec<u32>, BTreeSet<u32>)> = Vec::new();
    for (p_idx, reg) in regions {
        let footprint: &BTreeSet<u32> = match reg {
            Some(r) => r,
            None => {
                if colors.is_empty() {
                    colors.push((Vec::new(), BTreeSet::new()));
                }
                colors[0].0.push(*p_idx);
                continue;
            }
        };
        let k = colors
            .iter()
            .enumerate()
            .rev()
            .find(|(_, (_, occ))| !occ.is_disjoint(footprint))
            .map_or(0, |(j, _)| j + 1);
        if k == colors.len() {
            colors.push((Vec::new(), BTreeSet::new()));
        }
        colors[k].0.push(*p_idx);
        colors[k].1.extend(footprint.iter().copied());
    }
    colors.into_iter().map(|(pts, _)| pts).collect()
}

/// Build the Delaunay tetrahedralization by deterministic prefix
/// batches: conflict regions read-only across `threads` scoped
/// threads, application in EXACT BRIO order (batches never reorder).
/// The finished mesh is bitwise identical to the sequential kernel's
/// at any thread count — gated by the battery on general-position AND
/// degenerate inputs.
///
/// # Errors
/// Same surface as [`crate::delaunay::delaunay`].
///
/// # Panics
/// Only on kernel programmer contracts (a live real tet always
/// exists once bootstrapped).
pub fn delaunay_colored(
    points: &[fs_geom::Point3],
    threads: usize,
    window: usize,
    cx: &Cx<'_>,
) -> Result<(Tetrahedralization, ColoredStats), MeshError> {
    let (mut mesh, quad, order) = bootstrap_mesh(points)?;
    let work: Vec<u32> = order.iter().copied().filter(|i| !quad.contains(i)).collect();
    let threads = threads.max(1);
    let window = window.max(1);
    let mut stats = ColoredStats {
        threads: threads as u64,
        points: work.len() as u64,
        ..ColoredStats::default()
    };
    for win in work.chunks(window) {
        // Phase A: read-only conflict regions, deterministic thread
        // partition (contiguous chunks reassembled by position — the
        // schedule cannot change the result, only the wall clock).
        let chunk = win.len().div_ceil(threads);
        let mesh_ref = &mesh;
        let mut regions: Vec<(u32, Option<BTreeSet<u32>>)> = Vec::with_capacity(win.len());
        std::thread::scope(|scope| {
            let handles: Vec<_> = win
                .chunks(chunk)
                .map(|part| {
                    scope.spawn(move || {
                        part.iter()
                            .map(|&p| (p, conflict_region(mesh_ref, p)))
                            .collect::<Vec<_>>()
                    })
                })
                .collect();
            for h in handles {
                regions.extend(h.join().expect("region thread panicked"));
            }
        });
        // Phase B: flip-safe coloring.
        let colors = assign_colors(&regions);
        for c in &colors {
            stats.batches += 1;
            stats.largest_batch = stats.largest_batch.max(c.len() as u64);
            if c.len() == 1 {
                stats.singleton_batches += 1;
            }
        }
        // Phase C: color-by-color application, each color in BRIO
        // order (flips only across disjoint pairs).
        for color in &colors {
            for &p in color {
                mesh.insert(p);
            }
        }
        cx.checkpoint()?;
    }
    mesh.stats.tets_final = (0..mesh.tets.len() as u32)
        .filter(|&t| mesh.alive[t as usize] && !mesh.is_ghost(t))
        .count() as u64;
    let steiner_from = u32::try_from(points.len()).expect("point count fits u32");
    Ok((
        Tetrahedralization {
            mesh,
            steiner_from,
        },
        stats,
    ))
}

/// Apply each prefix batch REVERSED — the adversarial commutativity
/// probe: if batch regions were not truly pairwise disjoint, the
/// reversed order would produce a different mesh (compared
/// canonically, since allocation order legitimately differs).
///
/// # Errors
/// Same surface as [`delaunay_colored`].
pub fn delaunay_colored_reversed(
    points: &[fs_geom::Point3],
    window: usize,
    cx: &Cx<'_>,
) -> Result<Tetrahedralization, MeshError> {
    let (mut mesh, quad, order) = bootstrap_mesh(points)?;
    let work: Vec<u32> = order.iter().copied().filter(|i| !quad.contains(i)).collect();
    let window = window.max(1);
    for win in work.chunks(window.max(1)) {
        let mesh_ref = &mesh;
        let regions: Vec<(u32, Option<BTreeSet<u32>>)> = win
            .iter()
            .map(|&p| (p, conflict_region(mesh_ref, p)))
            .collect();
        let colors = assign_colors(&regions);
        for color in &colors {
            for &p in color.iter().rev() {
                mesh.insert(p);
            }
        }
        cx.checkpoint()?;
    }
    mesh.stats.tets_final = (0..mesh.tets.len() as u32)
        .filter(|&t| mesh.alive[t as usize] && !mesh.is_ghost(t))
        .count() as u64;
    let steiner_from = u32::try_from(points.len()).expect("point count fits u32");
    Ok(Tetrahedralization {
        mesh,
        steiner_from,
    })
}

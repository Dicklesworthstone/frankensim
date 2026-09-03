//! The polygon-soup repair suite (plan §7.2): duplicate/degenerate-face
//! removal, orientation unification (flood fill + signed-volume vote),
//! and small-hole fan filling — every action logged as a structured
//! RECEIPT (defect type, location, action) in the format the fs-io
//! quarantine consumes. Self-intersection FLAGGING is the validity-
//! certificates bead's (it needs fs-ivl broad/narrow phases to be
//! certified); this suite repairs what winding-robust queries cannot
//! absorb.

use crate::winding::Soup;
use std::collections::BTreeMap;
use std::fmt::Write as _;

/// One repair action's receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct RepairReceipt {
    /// Defect class ("duplicate-face", "degenerate-face",
    /// "flipped-patch", "boundary-hole").
    pub defect: &'static str,
    /// Where (face index or vertex list, rendered).
    pub location: String,
    /// What was done ("removed", "flipped n faces", "fan-filled").
    pub action: String,
}

impl RepairReceipt {
    /// Canonical JSON object.
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"defect\":\"{}\",\"location\":\"{}\",\"action\":\"{}\"}}",
            self.defect, self.location, self.action
        )
    }
}

/// The repair outcome: the healed soup + receipts.
#[derive(Debug)]
pub struct RepairOutcome {
    /// The repaired soup.
    pub soup: Soup,
    /// Per-defect receipts (deterministic order).
    pub receipts: Vec<RepairReceipt>,
}

impl RepairOutcome {
    /// The receipts as one JSON array (the quarantine log line).
    #[must_use]
    pub fn receipts_json(&self) -> String {
        let mut s = String::from("[");
        for (i, r) in self.receipts.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            let _ = write!(s, "{}", r.to_json());
        }
        s.push(']');
        s
    }
}

/// Run the full pipeline: dedupe → drop degenerates → unify orientation →
/// fan-fill boundary holes up to `max_hole_edges`.
#[must_use]
pub fn repair(mut soup: Soup, max_hole_edges: usize) -> RepairOutcome {
    let mut receipts = Vec::new();
    dedupe_faces(&mut soup, &mut receipts);
    drop_degenerates(&mut soup, &mut receipts);
    unify_orientation(&mut soup, &mut receipts);
    fill_holes(&mut soup, max_hole_edges, &mut receipts);
    RepairOutcome { soup, receipts }
}

fn canonical(tri: [u32; 3]) -> [u32; 3] {
    // Rotation-canonical, orientation-PRESERVING key start; combined with
    // the sorted key for duplicate detection regardless of orientation.
    let mut t = tri;
    t.sort_unstable();
    t
}

fn dedupe_faces(soup: &mut Soup, receipts: &mut Vec<RepairReceipt>) {
    let mut seen: BTreeMap<[u32; 3], usize> = BTreeMap::new();
    let mut keep = Vec::with_capacity(soup.triangles.len());
    for (i, &tri) in soup.triangles.iter().enumerate() {
        let key = canonical(tri);
        if let Some(&first) = seen.get(&key) {
            receipts.push(RepairReceipt {
                defect: "duplicate-face",
                location: format!("face {i} duplicates face {first} (vertices {key:?})"),
                action: "removed".to_string(),
            });
        } else {
            seen.insert(key, i);
            keep.push(tri);
        }
    }
    soup.triangles = keep;
}

fn drop_degenerates(soup: &mut Soup, receipts: &mut Vec<RepairReceipt>) {
    let positions = soup.positions.clone();
    let mut keep = Vec::with_capacity(soup.triangles.len());
    for (i, &tri) in soup.triangles.iter().enumerate() {
        let [a, b, c] = tri.map(|v| positions[v as usize]);
        let cross = cross3(b.delta_from(a), c.delta_from(a));
        let repeated = tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2];
        if repeated || cross.norm() < 1e-30 {
            receipts.push(RepairReceipt {
                defect: "degenerate-face",
                location: format!("face {i} (vertices {tri:?})"),
                action: "removed".to_string(),
            });
        } else {
            keep.push(tri);
        }
    }
    soup.triangles = keep;
}

fn cross3(a: fs_geom::Vec3, b: fs_geom::Vec3) -> fs_geom::Vec3 {
    fs_geom::Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}

/// Flood-fill orientation unification: neighbors sharing an edge must
/// traverse it in OPPOSITE directions; the shell is then globally oriented
/// outward by the sign of its enclosed volume.
fn unify_orientation(soup: &mut Soup, receipts: &mut Vec<RepairReceipt>) {
    let nt = soup.triangles.len();
    // Edge → incident (face, direction) map.
    let mut edge_faces: BTreeMap<[u32; 2], Vec<(usize, bool)>> = BTreeMap::new();
    for (fi, tri) in soup.triangles.iter().enumerate() {
        for c in 0..3 {
            let (a, b) = (tri[c], tri[(c + 1) % 3]);
            let key = if a < b { [a, b] } else { [b, a] };
            edge_faces.entry(key).or_default().push((fi, a < b));
        }
    }
    let mut flipped = vec![false; nt];
    let mut visited = vec![false; nt];
    let mut total_flips = 0usize;
    for seed in 0..nt {
        if visited[seed] {
            continue;
        }
        let mut stack = vec![seed];
        visited[seed] = true;
        while let Some(f) = stack.pop() {
            let tri = soup.triangles[f];
            for c in 0..3 {
                let (a, b) = (tri[c], tri[(c + 1) % 3]);
                let key = if a < b { [a, b] } else { [b, a] };
                let dir_here = (a < b) != flipped[f];
                for &(g, dir_there_raw) in &edge_faces[&key] {
                    if g == f || visited[g] {
                        continue;
                    }
                    // Consistent orientation = opposite traversal.
                    let needs_flip = (dir_there_raw != flipped[g]) == dir_here;
                    flipped[g] = needs_flip;
                    if needs_flip {
                        total_flips += 1;
                    }
                    visited[g] = true;
                    stack.push(g);
                }
            }
        }
    }
    for (f, &flip) in flipped.iter().enumerate() {
        if flip {
            soup.triangles[f].swap(1, 2);
        }
    }
    // Global orientation: the SIGNED enclosed volume of the now-consistent
    // shell (divergence theorem), negative iff every face points inward.
    // It is +-V for any closed shell, so the decision is never made from
    // rounding. The vertex-centroid winding number used before is only
    // meaningful when that centroid lies inside the solid; for a finned
    // heatsink, a ring or an L it lies in the air, the winding is a
    // floating-point residual of order 1e-16 with an arbitrary sign, and the
    // rotated four-fin heatsink (MEASURED 2026-09-02, fs-io's f32 weld) was
    // inverted by it — every downstream consumer then saw an inward surface
    // (fs-mesh's audit refused: retained volume +V against a closed-surface
    // volume of -V).
    let signed_volume: f64 = soup
        .triangles
        .iter()
        .map(|&[a, b, c]| {
            let (a, b, c) = (
                soup.positions[a as usize],
                soup.positions[b as usize],
                soup.positions[c as usize],
            );
            let bc = [
                b.y * c.z - b.z * c.y,
                b.z * c.x - b.x * c.z,
                b.x * c.y - b.y * c.x,
            ];
            (a.x * bc[0] + a.y * bc[1] + a.z * bc[2]) / 6.0
        })
        .sum();
    if std::env::var_os("FS_REPMESH_TRACE_ORIENT").is_some() {
        // Diagnostic: the old vote next to the new one, per shell.
        let n = soup.positions.len().max(1) as f64;
        let mut c = [0.0f64; 3];
        for p in &soup.positions {
            c = [c[0] + p.x / n, c[1] + p.y / n, c[2] + p.z / n];
        }
        let centroid_winding =
            crate::winding::winding_exact(soup, fs_geom::Point3::new(c[0], c[1], c[2]));
        eprintln!(
            "TRACE orient: {} faces, {} vertices: signed volume {signed_volume:e}, vertex-centroid winding {centroid_winding:e} (old rule would {})",
            soup.triangles.len(),
            soup.positions.len(),
            if centroid_winding < 0.0 {
                "flip"
            } else {
                "keep"
            }
        );
    }
    if signed_volume < 0.0 {
        for tri in &mut soup.triangles {
            tri.swap(1, 2);
        }
        total_flips += soup.triangles.len();
        receipts.push(RepairReceipt {
            defect: "flipped-patch",
            location: "global".to_string(),
            action: format!("inverted every face (signed enclosed volume was {signed_volume:e})"),
        });
    }
    if total_flips > 0 {
        receipts.push(RepairReceipt {
            defect: "flipped-patch",
            location: format!("{total_flips} faces"),
            action: "flipped to the component-consistent orientation".to_string(),
        });
    }
}

/// Fan-fill boundary loops of up to `max_hole_edges` edges.
fn fill_holes(soup: &mut Soup, max_hole_edges: usize, receipts: &mut Vec<RepairReceipt>) {
    // Boundary edges: traversed exactly once.
    let mut edge_count: BTreeMap<[u32; 2], (usize, (u32, u32))> = BTreeMap::new();
    for tri in &soup.triangles {
        for c in 0..3 {
            let (a, b) = (tri[c], tri[(c + 1) % 3]);
            let key = if a < b { [a, b] } else { [b, a] };
            let e = edge_count.entry(key).or_insert((0, (a, b)));
            e.0 += 1;
        }
    }
    // Directed boundary edges dest->... build loop maps (b -> a means the
    // hole is traversed a -> b on the missing side).
    let mut next_of: BTreeMap<u32, u32> = BTreeMap::new();
    for &(count, (a, b)) in edge_count.values() {
        if count == 1 {
            // The face traverses a->b; the hole's loop runs b->a.
            next_of.insert(b, a);
        }
    }
    let mut visited: std::collections::BTreeSet<u32> = std::collections::BTreeSet::new();
    let starts: Vec<u32> = next_of.keys().copied().collect();
    for start in starts {
        if visited.contains(&start) {
            continue;
        }
        // Walk the loop.
        let mut loop_verts = vec![start];
        visited.insert(start);
        let mut cur = start;
        loop {
            let Some(&nxt) = next_of.get(&cur) else {
                loop_verts.clear();
                break;
            };
            if nxt == start {
                break;
            }
            if !visited.insert(nxt) {
                loop_verts.clear();
                break;
            }
            loop_verts.push(nxt);
            cur = nxt;
        }
        if loop_verts.len() >= 3 && loop_verts.len() <= max_hole_edges {
            let anchor = loop_verts[0];
            for w in loop_verts[1..].windows(2) {
                soup.triangles.push([anchor, w[0], w[1]]);
            }
            receipts.push(RepairReceipt {
                defect: "boundary-hole",
                location: format!("loop of {} edges at vertex {anchor}", loop_verts.len()),
                action: format!("fan-filled with {} triangles", loop_verts.len() - 2),
            });
        }
    }
}

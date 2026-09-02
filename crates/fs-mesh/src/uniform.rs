//! Uniform 1→8 tetrahedral refinement (the "red" refinement of Bey 1995 /
//! Zhang 1995 / Liu–Joe 1996).
//!
//! Every tet is split by its six edge midpoints into four corner tets (each a
//! half-scale copy of the parent, listed in the parent's vertex order) and
//! four interior tets from the central octahedron cut along its SHORTEST
//! diagonal (ties broken by canonical vertex order). Liu and Joe (Math. Comp.
//! 65, 1996) prove that the shortest-diagonal choice keeps the quality of
//! every descendant bounded below by a constant times the parent's, so the
//! refinement can be applied repeatedly as an h-ladder without degenerating;
//! `tests/uniform_refine.rs` measures it. Interior children are oriented by
//! the exact predicate, never by a sign heuristic.
//!
//! Conforming by construction: a face shared by two parents is split by the
//! same three midpoints on both sides. Recovered walls (source faces) split
//! four ways in their own plane with the parent facet inherited, so boundary
//! classification by parent facet survives every rung. Region labels
//! replicate; volume is preserved up to the correctly rounded midpoints.
//! No Delaunay property is claimed or needed — this is the convergence-study
//! control, not a quality improver.

use fs_ivl::{Sign, orient3d};
use std::collections::BTreeMap;

/// One uniform refinement of a tet complex given as arrays. Children of
/// parent tet `i` occupy slots `8i..8i+8`, so per-parent labels replicate by
/// index.
pub(crate) struct UniformSplit {
    pub(crate) positions: Vec<[f64; 3]>,
    pub(crate) tets: Vec<[u32; 4]>,
    pub(crate) source_faces: Vec<([u32; 3], u32)>,
    /// Edge midpoints appended after the input vertices.
    pub(crate) midpoints_added: usize,
}

struct Midpoints {
    by_edge: BTreeMap<(u32, u32), u32>,
    positions: Vec<[f64; 3]>,
}

impl Midpoints {
    fn of(&mut self, a: u32, b: u32) -> u32 {
        let key = if a < b { (a, b) } else { (b, a) };
        if let Some(&m) = self.by_edge.get(&key) {
            return m;
        }
        let pa = self.positions[a as usize];
        let pb = self.positions[b as usize];
        let m = u32::try_from(self.positions.len()).expect("refined vertex count fits u32");
        self.positions.push([
            f64::midpoint(pa[0], pb[0]),
            f64::midpoint(pa[1], pb[1]),
            f64::midpoint(pa[2], pb[2]),
        ]);
        self.by_edge.insert(key, m);
        m
    }
}

fn squared_distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
}

/// Give `tet` the parent's exact orientation sign: swap the last two
/// vertices when the exact predicate disagrees with `want`. Convention-free
/// (whatever sign the parent carries, its children carry). A degenerate
/// parent (`Sign::Zero`) leaves children as listed; the caller's audit and
/// quality census own that refusal.
fn oriented(points: &[[f64; 3]], tet: [u32; 4], want: Sign) -> [u32; 4] {
    let p = |v: u32| points[v as usize];
    let sign = orient3d(p(tet[0]), p(tet[1]), p(tet[2]), p(tet[3]));
    if want != Sign::Zero && sign != Sign::Zero && sign != want {
        [tet[0], tet[1], tet[3], tet[2]]
    } else {
        tet
    }
}

/// Split every tet 1→8 and every source face 1→4.
pub(crate) fn split_uniform(
    positions: &[[f64; 3]],
    tets: &[[u32; 4]],
    source_faces: &[([u32; 3], u32)],
) -> UniformSplit {
    let input_vertices = positions.len();
    let mut mids = Midpoints {
        by_edge: BTreeMap::new(),
        positions: positions.to_vec(),
    };
    let mut out: Vec<[u32; 4]> = Vec::with_capacity(tets.len().saturating_mul(8));
    for &[x0, x1, x2, x3] in tets {
        let want = {
            let p = |v: u32| mids.positions[v as usize];
            orient3d(p(x0), p(x1), p(x2), p(x3))
        };
        let m01 = mids.of(x0, x1);
        let m02 = mids.of(x0, x2);
        let m03 = mids.of(x0, x3);
        let m12 = mids.of(x1, x2);
        let m13 = mids.of(x1, x3);
        let m23 = mids.of(x2, x3);
        // Corner tets: the parent scaled by one half about each vertex, in
        // the parent's own order — positively oriented whenever the parent is.
        out.push([x0, m01, m02, m03]);
        out.push([m01, x1, m12, m13]);
        out.push([m02, m12, x2, m23]);
        out.push([m03, m13, m23, x3]);
        // The octahedron's three diagonals join midpoints of opposite edges.
        // Shortest first; the canonical (index-ordered) diagonal breaks ties
        // so the split is a pure function of the input.
        let diagonals: [((u32, u32), [u32; 4]); 3] = [
            ((m01, m23), [m02, m03, m13, m12]),
            ((m02, m13), [m01, m03, m23, m12]),
            ((m03, m12), [m01, m02, m23, m13]),
        ];
        let mut best = 0usize;
        let mut best_len = f64::INFINITY;
        for (index, ((p, q), _)) in diagonals.iter().enumerate() {
            let len = squared_distance(mids.positions[*p as usize], mids.positions[*q as usize]);
            if len < best_len {
                best_len = len;
                best = index;
            }
        }
        let ((p, q), equator) = diagonals[best];
        for k in 0..4 {
            let e0 = equator[k];
            let e1 = equator[(k + 1) % 4];
            out.push(oriented(&mids.positions, [p, q, e0, e1], want));
        }
    }
    let mut faces = Vec::with_capacity(source_faces.len().saturating_mul(4));
    for &([a, b, c], parent) in source_faces {
        let ab = mids.of(a, b);
        let bc = mids.of(b, c);
        let ca = mids.of(c, a);
        faces.push(([a, ab, ca], parent));
        faces.push(([ab, b, bc], parent));
        faces.push(([ca, bc, c], parent));
        faces.push(([ab, bc, ca], parent));
    }
    let midpoints_added = mids.positions.len() - input_vertices;
    UniformSplit {
        positions: mids.positions,
        tets: out,
        source_faces: faces,
        midpoints_added,
    }
}

#[cfg(test)]
mod tests {
    use super::split_uniform;

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
    fn norm(a: [f64; 3]) -> f64 {
        dot(a, a).sqrt()
    }
    fn signed_volume(p: &[[f64; 3]; 4]) -> f64 {
        dot(sub(p[1], p[0]), cross(sub(p[2], p[0]), sub(p[3], p[0]))) / 6.0
    }
    fn min_dihedral_deg(p: &[[f64; 3]; 4]) -> f64 {
        let mut worst = 180.0f64;
        for (i, j) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] {
            let others: Vec<usize> = (0..4).filter(|&k| k != i && k != j).collect();
            let (k, l) = (others[0], others[1]);
            let n1 = cross(sub(p[j], p[i]), sub(p[k], p[i]));
            let n2 = cross(sub(p[l], p[i]), sub(p[j], p[i]));
            let cosine = -dot(n1, n2) / (norm(n1) * norm(n2));
            worst = worst.min(cosine.clamp(-1.0, 1.0).acos().to_degrees());
        }
        worst
    }
    fn census(positions: &[[f64; 3]], tets: &[[u32; 4]]) -> (f64, f64) {
        let mut volume = 0.0;
        let mut min_dih = 180.0f64;
        for t in tets {
            let p = [
                positions[t[0] as usize],
                positions[t[1] as usize],
                positions[t[2] as usize],
                positions[t[3] as usize],
            ];
            let v = signed_volume(&p);
            assert!(v > 0.0, "child is positively oriented, found {v}");
            volume += v;
            min_dih = min_dih.min(min_dihedral_deg(&p));
        }
        (volume, min_dih)
    }

    fn generations(start: [[f64; 3]; 4], rungs: usize) -> Vec<(usize, f64, f64)> {
        let mut positions = start.to_vec();
        let mut tets = vec![[0u32, 1, 2, 3]];
        let mut rows = Vec::new();
        for _ in 0..rungs {
            let split = split_uniform(&positions, &tets, &[]);
            positions = split.positions;
            tets = split.tets;
            let (volume, min_dih) = census(&positions, &tets);
            rows.push((tets.len(), volume, min_dih));
        }
        rows
    }

    #[test]
    fn regular_tet_children_conserve_volume_and_keep_dihedrals_for_five_generations() {
        let s = 0.5f64.sqrt();
        let start = [
            [1.0, 0.0, -s],
            [-1.0, 0.0, -s],
            [0.0, -1.0, s],
            [0.0, 1.0, s],
        ];
        let parent = signed_volume(&start);
        assert!(parent > 0.0);
        let rows = generations(start, 5);
        for (k, (count, volume, _)) in rows.iter().enumerate() {
            assert_eq!(*count, 8usize.pow(k as u32 + 1));
            assert!(
                ((volume - parent) / parent).abs() < 1e-11,
                "generation {k} volume {volume}"
            );
        }
        // Liu–Joe: quality bounded below for every generation. MEASURED
        // 2026-09-02: the regular tet's first generation drops from 70.53° to
        // the interior class's minimum, and no later generation goes lower.
        let first = rows[0].2;
        for (k, (_, _, min_dih)) in rows.iter().enumerate().skip(1) {
            assert!(
                *min_dih >= first - 1e-9,
                "generation {k} min dihedral {min_dih} below generation 1's {first}"
            );
        }
        assert!(
            first > 30.0,
            "regular tet children keep dihedrals above 30°, found {first}"
        );
    }

    #[test]
    fn a_skewed_tet_stays_bounded_over_five_generations() {
        let start = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.2, 0.7, 0.0],
            [0.3, 0.1, 0.35],
        ];
        let base = min_dihedral_deg(&start);
        let rows = generations(start, 5);
        let worst = rows.iter().map(|r| r.2).fold(180.0, f64::min);
        // Bounded, not monotone: later generations may not fall below the
        // worst class that has already appeared by the second generation.
        let by_second = rows[1].2;
        for (k, (_, _, min_dih)) in rows.iter().enumerate().skip(2) {
            assert!(
                *min_dih >= by_second - 1e-9,
                "generation {k} min dihedral {min_dih} below generation 2's {by_second}"
            );
        }
        assert!(worst >= 0.5 * base, "worst {worst} vs base {base}");
    }

    #[test]
    fn source_faces_split_four_ways_on_their_parent_and_stay_boundary_faces() {
        let start = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let faces = [([0u32, 1, 2], 7u32), ([0, 1, 3], 8)];
        let split = split_uniform(&start, &[[0, 1, 2, 3]], &faces);
        assert_eq!(split.tets.len(), 8);
        assert_eq!(split.source_faces.len(), 8);
        assert_eq!(split.midpoints_added, 6);
        let mut count = std::collections::BTreeMap::new();
        for t in &split.tets {
            for f in [
                [t[0], t[1], t[2]],
                [t[0], t[1], t[3]],
                [t[0], t[2], t[3]],
                [t[1], t[2], t[3]],
            ] {
                let mut key = f;
                key.sort_unstable();
                *count.entry(key).or_insert(0u32) += 1;
            }
        }
        for (face, parent) in &split.source_faces {
            let mut key = *face;
            key.sort_unstable();
            assert_eq!(
                count.get(&key),
                Some(&1),
                "refined wall {face:?} is a boundary face"
            );
            assert!(matches!(parent, 7 | 8));
            // Coplanar with the parent: z == 0 for parent 7, y == 0 for parent 8.
            for &v in face {
                let p = split.positions[v as usize];
                if *parent == 7 {
                    assert_eq!(p[2], 0.0);
                } else {
                    assert_eq!(p[1], 0.0);
                }
            }
        }
    }
}

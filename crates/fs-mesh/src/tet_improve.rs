//! Quality improvement of a tet mesh whose boundary is frozen (IRP M1.3-W4).
//!
//! Curving a boundary (snapping vertices onto an exact surface) can fold the
//! thin tets beside it, and uniform refinement keeps whatever slivers the
//! volumetricizer left. This pass repairs both without touching geometry
//! that defines the domain: frozen vertices (the boundary) never move and
//! boundary faces are never flipped, so the union, every boundary face and
//! the enclosed volume are unchanged.
//!
//! Two local operations, accepted only when they strictly raise the worst
//! quality of the tets they touch:
//! - smart smoothing of a free vertex over a few deterministic candidates
//!   (the centroid of its neighbours and the ideal apex over the opposite
//!   face of its worst tet, with step halving); an inverted star is first
//!   untangled by its normalized signed volume;
//! - 2-3 flips across interior faces and 3-2 flips around interior edges of
//!   degree three (Klingner and Shewchuk 2007 use the same moves with more
//!   aggressive schedules).
//!
//! Orientation is judged against the caller's sign convention (fs-mesh tets
//! carry one consistent sign) and re-checked exactly at the end; any tet left
//! inverted refuses. Quality is not certified: the result reports what it
//! reached.

use fs_exec::Cx;
use fs_ivl::{Sign, orient3d};
use std::collections::BTreeSet;

/// Improvement policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImproveOptions {
    /// Tets below this minimum dihedral angle [degrees] are worked on.
    pub target_min_dihedral_deg: f64,
    /// Sweeps over the work list.
    pub max_passes: u32,
}

/// What the pass did and reached.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ImproveStats {
    /// Inverted or zero-volume tets before.
    pub inverted_before: usize,
    /// Smallest minimum dihedral among valid tets before [degrees].
    pub min_dihedral_before_deg: f64,
    /// Smallest minimum dihedral after [degrees].
    pub min_dihedral_after_deg: f64,
    /// Tets still below the target after.
    pub below_target_after: usize,
    /// Accepted vertex moves.
    pub smoothed: usize,
    /// Accepted 2-3 flips.
    pub flips_23: usize,
    /// Accepted 3-2 flips.
    pub flips_32: usize,
    /// Sweeps run.
    pub passes: u32,
}

impl ImproveStats {
    /// Canonical JSON object.
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"inverted_before\":{},\"min_dihedral_before_deg\":{:e},\
             \"min_dihedral_after_deg\":{:e},\"below_target_after\":{},\"smoothed\":{},\
             \"flips_23\":{},\"flips_32\":{},\"passes\":{}}}",
            self.inverted_before,
            self.min_dihedral_before_deg,
            self.min_dihedral_after_deg,
            self.below_target_after,
            self.smoothed,
            self.flips_23,
            self.flips_32,
            self.passes
        )
    }
}

/// Refusals of [`improve_tets`].
#[derive(Debug, Clone, PartialEq)]
pub enum ImproveError {
    /// Options are out of range or the sign is not ±1.
    InvalidInput(&'static str),
    /// A tet is still inverted or flat after the pass.
    Inverted {
        /// Tet index in the output.
        tet: usize,
    },
    /// The oriented volume sum moved: the tets no longer tile the domain.
    VolumeNotConserved {
        /// Oriented volume before [length^3].
        before: f64,
        /// Oriented volume after [length^3].
        after: f64,
    },
    /// Cancellation was observed.
    Cancelled,
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

fn corners(positions: &[[f64; 3]], t: [u32; 4]) -> [[f64; 3]; 4] {
    t.map(|v| positions[v as usize])
}

/// Valid tets score their minimum dihedral in `[0, 70.6]` degrees; inverted
/// or flat ones score their normalized signed volume minus one, below every
/// valid score, so maximizing it untangles.
fn quality(p: [[f64; 3]; 4], sign: f64) -> f64 {
    let volume = sign * dot(sub(p[1], p[0]), cross(sub(p[2], p[0]), sub(p[3], p[0])));
    if volume > 0.0 {
        return crate::exude::min_dihedral_deg(&p);
    }
    let longest = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        .iter()
        .map(|&(a, b)| dot(sub(p[a], p[b]), sub(p[a], p[b])))
        .fold(0.0, f64::max);
    let scale = longest * fs_math::det::sqrt(longest);
    if scale > 0.0 {
        (volume / scale).max(-1.0) - 1.0
    } else {
        -2.0
    }
}

fn signed_volume(p: [[f64; 3]; 4]) -> f64 {
    dot(sub(p[1], p[0]), cross(sub(p[2], p[0]), sub(p[3], p[0]))) / 6.0
}

/// `t` reordered to start with `lead`, as an even permutation of `t` (the
/// last two slots swap when needed), so it keeps `t`'s orientation.
fn arrange(t: [u32; 4], lead: &[u32]) -> [u32; 4] {
    let mut order: Vec<usize> = lead
        .iter()
        .map(|v| t.iter().position(|w| w == v).expect("lead vertex is in t"))
        .collect();
    for i in 0..4 {
        if !order.contains(&i) {
            order.push(i);
        }
    }
    let inversions = (0..4)
        .flat_map(|i| (i + 1..4).map(move |j| (i, j)))
        .filter(|&(i, j)| order[i] > order[j])
        .count();
    let mut out = [t[order[0]], t[order[1]], t[order[2]], t[order[3]]];
    if inversions % 2 == 1 {
        out.swap(2, 3);
    }
    out
}

fn exact_valid(positions: &[[f64; 3]], t: [u32; 4], want: Sign) -> bool {
    let [a, b, c, d] = corners(positions, t);
    orient3d(a, b, c, d) == want
}

struct Mesh<'a> {
    positions: &'a mut [[f64; 3]],
    tets: Vec<[u32; 4]>,
    alive: Vec<bool>,
    incident: Vec<Vec<u32>>,
    frozen: &'a [bool],
    boundary_faces: BTreeSet<[u32; 3]>,
    sign: f64,
}

fn sorted3(mut f: [u32; 3]) -> [u32; 3] {
    f.sort_unstable();
    f
}

impl Mesh<'_> {
    fn q(&self, t: usize) -> f64 {
        quality(corners(self.positions, self.tets[t]), self.sign)
    }

    fn star_quality(&self, v: u32) -> f64 {
        self.incident[v as usize]
            .iter()
            .map(|&t| self.q(t as usize))
            .fold(f64::INFINITY, f64::min)
    }

    fn detach(&mut self, t: usize) {
        for v in self.tets[t] {
            self.incident[v as usize].retain(|&x| x as usize != t);
        }
    }

    fn attach(&mut self, t: usize) {
        for v in self.tets[t] {
            self.incident[v as usize].push(t as u32);
        }
    }

    /// Tets containing every vertex of `vs`, ascending.
    fn sharing(&self, vs: &[u32]) -> Vec<usize> {
        let mut out: Vec<usize> = self.incident[vs[0] as usize]
            .iter()
            .map(|&t| t as usize)
            .filter(|&t| vs.iter().all(|v| self.tets[t].contains(v)))
            .collect();
        out.sort_unstable();
        out
    }

    fn smooth(&mut self, v: u32) -> bool {
        if self.frozen[v as usize] || self.incident[v as usize].is_empty() {
            return false;
        }
        let star: Vec<usize> = self.incident[v as usize]
            .iter()
            .map(|&t| t as usize)
            .collect();
        let before = self.star_quality(v);
        let origin = self.positions[v as usize];
        let neighbours: BTreeSet<u32> = star
            .iter()
            .flat_map(|&t| self.tets[t])
            .filter(|&w| w != v)
            .collect();
        let n = neighbours.len() as f64;
        let centroid: [f64; 3] = core::array::from_fn(|k| {
            neighbours
                .iter()
                .map(|&w| self.positions[w as usize][k])
                .sum::<f64>()
                / n
        });
        // Ideal apex over the face opposite `v` in its worst tet.
        let worst = star
            .iter()
            .copied()
            .min_by(|&a, &b| self.q(a).total_cmp(&self.q(b)).then(a.cmp(&b)))
            .expect("non-empty star");
        let face: Vec<[f64; 3]> = self.tets[worst]
            .iter()
            .filter(|&&w| w != v)
            .map(|&w| self.positions[w as usize])
            .collect();
        let fc: [f64; 3] = core::array::from_fn(|k| (face[0][k] + face[1][k] + face[2][k]) / 3.0);
        let mut normal = cross(sub(face[1], face[0]), sub(face[2], face[0]));
        let len = fs_math::det::sqrt(dot(normal, normal));
        let edge = (fs_math::det::sqrt(dot(sub(face[1], face[0]), sub(face[1], face[0])))
            + fs_math::det::sqrt(dot(sub(face[2], face[1]), sub(face[2], face[1])))
            + fs_math::det::sqrt(dot(sub(face[0], face[2]), sub(face[0], face[2]))))
            / 3.0;
        let mut targets = vec![centroid];
        if len > 0.0 {
            normal = normal.map(|x| x / len);
            // The valid side is where the tet with `v` restored has the mesh sign.
            let mut probe = self.tets[worst];
            let slot = probe.iter().position(|&w| w == v).expect("v in its star");
            let height = (2.0f64 / 3.0).sqrt() * edge;
            let up: [f64; 3] = core::array::from_fn(|k| fc[k] + height * normal[k]);
            self.positions[v as usize] = up;
            probe[slot] = v;
            let up_ok = quality(corners(self.positions, probe), self.sign) > 0.0;
            self.positions[v as usize] = origin;
            let apex: [f64; 3] = if up_ok {
                up
            } else {
                core::array::from_fn(|k| fc[k] - height * normal[k])
            };
            targets.push(apex);
        }
        let mut best = (before, origin);
        for target in targets {
            let mut lambda = 1.0;
            for _ in 0..6 {
                let p: [f64; 3] =
                    core::array::from_fn(|k| origin[k] + lambda * (target[k] - origin[k]));
                self.positions[v as usize] = p;
                let score = self.star_quality(v);
                if score > best.0 {
                    best = (score, p);
                }
                lambda *= 0.5;
            }
        }
        self.positions[v as usize] = best.1;
        best.1 != origin
    }

    fn is_boundary_face(&self, f: [u32; 3]) -> bool {
        self.boundary_faces.contains(&sorted3(f))
    }

    /// 2-3 flip of `t` across its face opposite local vertex `k`. With the
    /// face ordered so `(f0, f1, f2, d)` keeps `t`'s orientation and `e` the
    /// neighbour's apex, the new tets `(f_i, f_i+1, e, d)` carry the mesh
    /// orientation exactly when the flip is valid (segment `de` through the
    /// face), so the sign test below is the validity test.
    fn flip23(&mut self, t: usize, k: usize) -> bool {
        let tet = self.tets[t];
        let d = tet[k];
        let lead = arrange(tet, &[d]);
        // (d, x, y, z) is even, so (x, z, y, d) is even too.
        let f = [lead[1], lead[3], lead[2]];
        if self.is_boundary_face(f) {
            return false;
        }
        let others: Vec<usize> = self.sharing(&f).into_iter().filter(|&u| u != t).collect();
        let &[u] = others.as_slice() else {
            return false;
        };
        let e = *self.tets[u]
            .iter()
            .find(|w| !f.contains(w))
            .expect("neighbour apex");
        // The new edge d-e must not already exist.
        if !self.sharing(&[d, e]).is_empty() {
            return false;
        }
        let candidates = [[f[0], f[1], e, d], [f[1], f[2], e, d], [f[2], f[0], e, d]];
        let before = self.q(t).min(self.q(u));
        let after = candidates
            .iter()
            .map(|c| quality(corners(self.positions, *c), self.sign))
            .fold(f64::INFINITY, f64::min);
        if !(after > before && after > 0.0) {
            return false;
        }
        self.detach(t);
        self.detach(u);
        self.tets[t] = candidates[0];
        self.tets[u] = candidates[1];
        self.tets.push(candidates[2]);
        self.alive.push(true);
        self.attach(t);
        self.attach(u);
        self.attach(self.tets.len() - 1);
        true
    }

    /// 3-2 flip removing edge (a, b) when its ring has three tets. Each ring
    /// tet ordered as an even `(a, b, p, q)` gives the ring cycle p -> q; for
    /// the cycle P0 -> P1 -> P2 the new tets `(P0, P1, P2, b)` and
    /// `(P0, P2, P1, a)` carry the mesh orientation exactly when the flip is
    /// valid (edge ab through the new face).
    fn flip32(&mut self, a: u32, b: u32) -> bool {
        let ring = self.sharing(&[a, b]);
        if ring.len() != 3 {
            return false;
        }
        let mut next = std::collections::BTreeMap::new();
        for &r in &ring {
            let [_, _, p, q] = arrange(self.tets[r], &[a, b]);
            // Every ring face through (a, b) must be interior.
            if self.is_boundary_face([a, b, p]) || self.is_boundary_face([a, b, q]) {
                return false;
            }
            if next.insert(p, q).is_some() {
                return false;
            }
        }
        let p0 = *next.keys().next().expect("three ring tets");
        let p1 = next[&p0];
        let Some(&p2) = next.get(&p1) else {
            return false;
        };
        if next.get(&p2) != Some(&p0) {
            return false;
        }
        // The new face must not already exist.
        if !self.sharing(&[p0, p1, p2]).is_empty() {
            return false;
        }
        let candidates = [[p0, p1, p2, b], [p0, p2, p1, a]];
        let before = ring
            .iter()
            .map(|&r| self.q(r))
            .fold(f64::INFINITY, f64::min);
        let after = candidates
            .iter()
            .map(|c| quality(corners(self.positions, *c), self.sign))
            .fold(f64::INFINITY, f64::min);
        if !(after > before && after > 0.0) {
            return false;
        }
        for &r in &ring {
            self.detach(r);
        }
        self.tets[ring[0]] = candidates[0];
        self.tets[ring[1]] = candidates[1];
        self.alive[ring[2]] = false;
        self.tets[ring[2]] = [u32::MAX; 4];
        self.attach(ring[0]);
        self.attach(ring[1]);
        true
    }
}

/// Improve `tets` in place, moving only vertices whose `frozen` flag is false
/// and flipping only across faces not in `boundary_faces`. `sign` is the
/// sign of `dot(b-a, (c-a) x (d-a))` every valid tet of this mesh carries.
///
/// # Errors
/// Invalid options, a tet still inverted or flat afterwards (exact
/// predicate), cancellation.
pub fn improve_tets(
    positions: &mut [[f64; 3]],
    tets: &mut Vec<[u32; 4]>,
    frozen: &[bool],
    boundary_faces: &[[u32; 3]],
    sign: f64,
    options: ImproveOptions,
    cx: &Cx<'_>,
) -> Result<ImproveStats, ImproveError> {
    if !(sign == 1.0 || sign == -1.0) {
        return Err(ImproveError::InvalidInput("sign must be +1 or -1"));
    }
    if !(options.target_min_dihedral_deg.is_finite()
        && (0.0..=60.0).contains(&options.target_min_dihedral_deg))
    {
        return Err(ImproveError::InvalidInput(
            "target dihedral must be in [0, 60] degrees",
        ));
    }
    if frozen.len() != positions.len() {
        return Err(ImproveError::InvalidInput("one frozen flag per vertex"));
    }
    let mut incident = vec![Vec::new(); positions.len()];
    for (i, t) in tets.iter().enumerate() {
        for &v in t {
            incident
                .get_mut(v as usize)
                .ok_or(ImproveError::InvalidInput("tet vertex out of range"))?
                .push(i as u32);
        }
    }
    let alive = vec![true; tets.len()];
    let oriented_volume = |positions: &[[f64; 3]], tets: &[[u32; 4]]| {
        tets.iter()
            .map(|t| sign * signed_volume(corners(positions, *t)))
            .sum::<f64>()
    };
    let volume_before = oriented_volume(positions, tets);
    let mut mesh = Mesh {
        positions,
        tets: std::mem::take(tets),
        alive,
        incident,
        frozen,
        boundary_faces: boundary_faces.iter().map(|f| sorted3(*f)).collect(),
        sign,
    };
    let mut stats = ImproveStats {
        min_dihedral_before_deg: 180.0,
        ..ImproveStats::default()
    };
    for t in 0..mesh.tets.len() {
        let q = mesh.q(t);
        if q > 0.0 {
            stats.min_dihedral_before_deg = stats.min_dihedral_before_deg.min(q);
        } else {
            stats.inverted_before += 1;
        }
    }
    let target = options.target_min_dihedral_deg;
    for pass in 0..options.max_passes {
        cx.checkpoint().map_err(|_| ImproveError::Cancelled)?;
        let bad: Vec<usize> = (0..mesh.tets.len())
            .filter(|&t| mesh.alive[t] && mesh.q(t) < target)
            .collect();
        if bad.is_empty() {
            break;
        }
        stats.passes = pass + 1;
        let mut changed = false;
        let vertices: BTreeSet<u32> = bad.iter().flat_map(|&t| mesh.tets[t]).collect();
        for (i, v) in vertices.into_iter().enumerate() {
            if i % 256 == 0 {
                cx.checkpoint().map_err(|_| ImproveError::Cancelled)?;
            }
            if mesh.smooth(v) {
                stats.smoothed += 1;
                changed = true;
            }
        }
        for t in bad {
            if !mesh.alive[t] || mesh.q(t) >= target {
                continue;
            }
            let mut done = false;
            for k in 0..4 {
                if mesh.flip23(t, k) {
                    stats.flips_23 += 1;
                    done = true;
                    break;
                }
            }
            if !done {
                let tet = mesh.tets[t];
                for (a, b) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] {
                    if mesh.flip32(tet[a], tet[b]) {
                        stats.flips_32 += 1;
                        done = true;
                        break;
                    }
                }
            }
            changed |= done;
        }
        if !changed {
            break;
        }
    }
    let Mesh {
        positions,
        tets: all,
        alive,
        ..
    } = mesh;
    let want = all
        .iter()
        .zip(&alive)
        .find(|(_, a)| **a)
        .map(|(t, _)| {
            let [a, b, c, d] = corners(positions, *t);
            orient3d(a, b, c, d)
        })
        .unwrap_or(Sign::Zero);
    *tets = all
        .into_iter()
        .zip(alive)
        .filter(|(_, a)| *a)
        .map(|(t, _)| t)
        .collect();
    let volume_after = oriented_volume(positions, tets);
    if !((volume_after - volume_before).abs() <= 1e-10 * volume_before.abs()) {
        return Err(ImproveError::VolumeNotConserved {
            before: volume_before,
            after: volume_after,
        });
    }
    stats.min_dihedral_after_deg = 180.0;
    for (i, t) in tets.iter().enumerate() {
        let q = quality(corners(positions, *t), sign);
        if !(q > 0.0) || !exact_valid(positions, *t, want) || want == Sign::Zero {
            return Err(ImproveError::Inverted { tet: i });
        }
        stats.min_dihedral_after_deg = stats.min_dihedral_after_deg.min(q);
        if q < target {
            stats.below_target_after += 1;
        }
    }
    Ok(stats)
}

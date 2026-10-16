//! Generalized winding numbers (plan §7.2, Barill et al.): ROBUST
//! inside/outside classification on broken input — "the single most
//! effective modern trick in mesh robustness". Exact evaluation sums
//! per-triangle solid angles (van Oosterom–Strackee); the octree DIPOLE
//! approximation aggregates far triangles into area-weighted normal
//! moments so million-triangle soups classify fast, with the
//! approximation error MEASURED against exact on fixtures (rmesh-004).

use fs_geom::{Point3, Vec3};

/// A triangle soup (no connectivity assumed — that is the point).
#[derive(Debug, Clone)]
pub struct Soup {
    /// Vertex positions.
    pub positions: Vec<Point3>,
    /// Triangles as vertex index triples.
    pub triangles: Vec<[u32; 3]>,
}

impl Soup {
    /// Triangle corner positions.
    #[must_use]
    pub fn tri(&self, t: usize) -> [Point3; 3] {
        let [a, b, c] = self.triangles[t];
        [
            self.positions[a as usize],
            self.positions[b as usize],
            self.positions[c as usize],
        ]
    }
}

/// Solid angle of triangle `[a, b, c]` seen from `p`, divided by 4π
/// (van Oosterom–Strackee — exact up to floating point).
#[must_use]
pub fn triangle_winding(p: Point3, a: Point3, b: Point3, c: Point3) -> f64 {
    let (va, vb, vc) = (a.delta_from(p), b.delta_from(p), c.delta_from(p));
    let (la, lb, lc) = (va.norm(), vb.norm(), vc.norm());
    let det = va.dot(cross(vb, vc));
    let denom = la * lb * lc + va.dot(vb) * lc + vb.dot(vc) * la + vc.dot(va) * lb;
    let omega = 2.0 * det.atan2(denom);
    omega / (4.0 * core::f64::consts::PI)
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}

/// Exact generalized winding number of the whole soup at `p`.
#[must_use]
pub fn winding_exact(soup: &Soup, p: Point3) -> f64 {
    (0..soup.triangles.len())
        .map(|t| {
            let [a, b, c] = soup.tri(t);
            triangle_winding(p, a, b, c)
        })
        .sum()
}

/// One octree node of the dipole hierarchy.
struct DipoleNode {
    /// Aggregate area-weighted normal (the dipole moment).
    moment: Vec3,
    /// Area-weighted centroid.
    centroid: Point3,
    /// Node bounding radius around the centroid.
    radius: f64,
    /// Triangles owned directly (leaves).
    tris: Vec<u32>,
    /// Child nodes (empty at leaves).
    children: Vec<DipoleNode>,
}

/// The dipole-accelerated winding evaluator (Barill et al. §3): far nodes
/// contribute the dipole term; near nodes recurse; leaves are exact. The
/// octree stores triangle INDICES only — callers pass the soup at
/// evaluation (no self-referential borrows; one octree, many queries).
pub struct WindingOctree {
    root: DipoleNode,
    /// Accuracy knob β: nodes farther than `β·radius` use the dipole term.
    beta: f64,
}

impl WindingOctree {
    /// Build over a soup with accuracy parameter `beta` (≥ 1; larger =
    /// more accurate, slower; 2.0 is the paper's sweet spot).
    #[must_use]
    pub fn build(soup: &Soup, beta: f64) -> Self {
        let tris: Vec<u32> = (0..soup.triangles.len() as u32).collect();
        let root = Self::build_node(soup, tris, 0);
        WindingOctree {
            root,
            beta: beta.max(1.0),
        }
    }

    fn tri_centroid_area_normal(soup: &Soup, t: u32) -> (Point3, f64, Vec3) {
        let [a, b, c] = soup.tri(t as usize);
        let n2 = cross(b.delta_from(a), c.delta_from(a));
        let area = 0.5 * n2.norm();
        let centroid = Point3::new(
            (a.x + b.x + c.x) / 3.0,
            (a.y + b.y + c.y) / 3.0,
            (a.z + b.z + c.z) / 3.0,
        );
        (centroid, area, n2.scale(0.5))
    }

    fn build_node(soup: &Soup, tris: Vec<u32>, depth: u32) -> DipoleNode {
        // Aggregate dipole data.
        let mut moment = Vec3::new(0.0, 0.0, 0.0);
        let mut wsum = 0.0f64;
        let mut csum = Vec3::new(0.0, 0.0, 0.0);
        for &t in &tris {
            let (c, area, n) = Self::tri_centroid_area_normal(soup, t);
            moment = Vec3::new(moment.x + n.x, moment.y + n.y, moment.z + n.z);
            let w = area.max(1e-300);
            wsum += w;
            csum = Vec3::new(
                c.x.mul_add(w, csum.x),
                c.y.mul_add(w, csum.y),
                c.z.mul_add(w, csum.z),
            );
        }
        let centroid = Point3::new(csum.x / wsum, csum.y / wsum, csum.z / wsum);
        let mut radius = 0.0f64;
        for &t in &tris {
            let corners = soup.tri(t as usize);
            for p in corners {
                radius = radius.max(p.delta_from(centroid).norm());
            }
        }
        if tris.len() <= 16 || depth >= 20 {
            return DipoleNode {
                moment,
                centroid,
                radius,
                tris,
                children: Vec::new(),
            };
        }
        // Octant split around the centroid (by triangle centroid).
        let mut buckets: [Vec<u32>; 8] = Default::default();
        for &t in &tris {
            let (c, _, _) = Self::tri_centroid_area_normal(soup, t);
            let idx = usize::from(c.x >= centroid.x)
                | (usize::from(c.y >= centroid.y) << 1)
                | (usize::from(c.z >= centroid.z) << 2);
            buckets[idx].push(t);
        }
        // Degenerate split (everything in one bucket): stay a leaf.
        if buckets.iter().filter(|b| !b.is_empty()).count() <= 1 {
            return DipoleNode {
                moment,
                centroid,
                radius,
                tris,
                children: Vec::new(),
            };
        }
        let children = buckets
            .into_iter()
            .filter(|b| !b.is_empty())
            .map(|b| Self::build_node(soup, b, depth + 1))
            .collect();
        DipoleNode {
            moment,
            centroid,
            radius,
            tris: Vec::new(),
            children,
        }
    }

    fn eval_node(&self, soup: &Soup, node: &DipoleNode, p: Point3) -> f64 {
        // Flux convention: w(p) = (1/4π) ∮ (x − p)·n dA / |x − p|³, so the
        // dipole term points FROM p TO the mass centroid.
        let d = node.centroid.delta_from(p);
        let dist = d.norm();
        if dist > self.beta * node.radius && dist > 0.0 {
            return node.moment.dot(d) / (4.0 * core::f64::consts::PI * dist * dist * dist);
        }
        if node.children.is_empty() {
            return node
                .tris
                .iter()
                .map(|&t| {
                    let [a, b, c] = soup.tri(t as usize);
                    triangle_winding(p, a, b, c)
                })
                .sum();
        }
        node.children
            .iter()
            .map(|ch| self.eval_node(soup, ch, p))
            .sum()
    }

    /// Approximate generalized winding number at `p`.
    #[must_use]
    pub fn winding(&self, soup: &Soup, p: Point3) -> f64 {
        self.eval_node(soup, &self.root, p)
    }

    /// Classify by the robust `w > 0.5` rule.
    #[must_use]
    pub fn inside(&self, soup: &Soup, p: Point3) -> bool {
        self.winding(soup, p) > 0.5
    }
}

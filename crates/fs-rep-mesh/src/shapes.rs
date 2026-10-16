//! Mesh fixture generators (PUBLIC — the mesh-side test vocabulary,
//! sibling to fs-geom's analytic fixtures): a unit cube soup, a
//! subdivided icosphere, and corruption helpers for the repair battery's
//! nightmare corpus.

use crate::winding::Soup;
use fs_geom::{Point3, Vec3};

/// Axis-aligned cube soup: 8 vertices, 12 outward-oriented triangles.
#[must_use]
pub fn cube(center: Point3, half: f64) -> Soup {
    let (c, h) = (center, half);
    let positions = vec![
        Point3::new(c.x - h, c.y - h, c.z - h), // 0
        Point3::new(c.x + h, c.y - h, c.z - h), // 1
        Point3::new(c.x + h, c.y + h, c.z - h), // 2
        Point3::new(c.x - h, c.y + h, c.z - h), // 3
        Point3::new(c.x - h, c.y - h, c.z + h), // 4
        Point3::new(c.x + h, c.y - h, c.z + h), // 5
        Point3::new(c.x + h, c.y + h, c.z + h), // 6
        Point3::new(c.x - h, c.y + h, c.z + h), // 7
    ];
    let triangles = vec![
        // -z (outward normal -z: CW seen from +z = CCW from -z)
        [0, 2, 1],
        [0, 3, 2],
        // +z
        [4, 5, 6],
        [4, 6, 7],
        // -y
        [0, 1, 5],
        [0, 5, 4],
        // +y
        [3, 7, 6],
        [3, 6, 2],
        // -x
        [0, 4, 7],
        [0, 7, 3],
        // +x
        [1, 2, 6],
        [1, 6, 5],
    ];
    Soup {
        positions,
        triangles,
    }
}

/// Icosphere: an icosahedron subdivided `subdivisions` times, vertices
/// projected to radius `radius` around `center`. Outward-oriented, closed.
#[must_use]
pub fn icosphere(center: Point3, radius: f64, subdivisions: u32) -> Soup {
    let phi = f64::midpoint(1.0, 5.0f64.sqrt());
    let base = [
        (-1.0, phi, 0.0),
        (1.0, phi, 0.0),
        (-1.0, -phi, 0.0),
        (1.0, -phi, 0.0),
        (0.0, -1.0, phi),
        (0.0, 1.0, phi),
        (0.0, -1.0, -phi),
        (0.0, 1.0, -phi),
        (phi, 0.0, -1.0),
        (phi, 0.0, 1.0),
        (-phi, 0.0, -1.0),
        (-phi, 0.0, 1.0),
    ];
    let mut positions: Vec<Point3> = base
        .iter()
        .map(|&(x, y, z)| project(center.offset(Vec3::new(x, y, z)), center, radius))
        .collect();
    let mut triangles: Vec<[u32; 3]> = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    for _ in 0..subdivisions {
        let mut midpoint: std::collections::BTreeMap<(u32, u32), u32> =
            std::collections::BTreeMap::new();
        let mut next = Vec::with_capacity(triangles.len() * 4);
        for &[a, b, c] in &triangles {
            let mut mid = |u: u32, v: u32, positions: &mut Vec<Point3>| -> u32 {
                let key = (u.min(v), u.max(v));
                if let Some(&m) = midpoint.get(&key) {
                    return m;
                }
                let (pu, pv) = (positions[u as usize], positions[v as usize]);
                let m = project(
                    Point3::new(
                        f64::midpoint(pu.x, pv.x),
                        f64::midpoint(pu.y, pv.y),
                        f64::midpoint(pu.z, pv.z),
                    ),
                    center,
                    radius,
                );
                positions.push(m);
                let idx = (positions.len() - 1) as u32;
                midpoint.insert(key, idx);
                idx
            };
            let ab = mid(a, b, &mut positions);
            let bc = mid(b, c, &mut positions);
            let ca = mid(c, a, &mut positions);
            next.extend_from_slice(&[[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]);
        }
        triangles = next;
    }
    Soup {
        positions,
        triangles,
    }
}

fn project(p: Point3, center: Point3, radius: f64) -> Point3 {
    // Direction from the CENTER (a latent off-origin bug — projecting
    // the absolute vector — was caught by fs-topo's exact
    // self-intersection certificate: off-center subdivided icospheres
    // came out spiky and genuinely self-intersecting).
    let d = p.delta_from(center);
    let n = d.norm().max(1e-300);
    center.offset(d.scale(radius / n))
}

/// Corrupt a soup for the repair battery: duplicate `dups` faces, insert
/// `degens` degenerate faces, flip the orientation of every face in
/// `flip_range`, and punch a hole by deleting face `punch`. Deterministic.
#[must_use]
pub fn corrupt(
    mut soup: Soup,
    dups: usize,
    degens: usize,
    flip_range: core::ops::Range<usize>,
    punch: Option<usize>,
) -> Soup {
    for i in 0..dups.min(soup.triangles.len()) {
        let t = soup.triangles[i * 7 % soup.triangles.len()];
        soup.triangles.push(t);
    }
    for i in 0..degens {
        let v = (i % soup.positions.len()) as u32;
        soup.triangles
            .push([v, v, (v + 1) % soup.positions.len() as u32]);
    }
    for f in flip_range {
        if f < soup.triangles.len() {
            soup.triangles[f].swap(1, 2);
        }
    }
    if let Some(p) = punch
        && p < soup.triangles.len()
    {
        soup.triangles.remove(p);
    }
    soup
}

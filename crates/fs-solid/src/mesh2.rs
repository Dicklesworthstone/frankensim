//! Structured 2D body-fitted meshes: P1 triangles and Q1 quads on
//! rectangles, and mapped Q1 panels for benchmark geometries (Cook's
//! membrane). Deliberately minimal — the unstructured 3D pipeline is
//! fs-mesh's bead; this module exists so the elasticity weak forms
//! have a body-fitted frontend with named boundary patches.

/// A named boundary patch (rectangles and mapped panels tag their four
/// sides).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Patch {
    /// x = min side (mapped: the left edge).
    Left,
    /// x = max side.
    Right,
    /// y = min side.
    Bottom,
    /// y = max side.
    Top,
}

/// A 2D mesh of P1 triangles or Q1 quads (one element type per mesh).
#[derive(Debug, Clone)]
pub struct Mesh2 {
    /// Node coordinates.
    pub nodes: Vec<[f64; 2]>,
    /// Element connectivity: 3 nodes (triangles) or 4 nodes (quads,
    /// counterclockwise).
    pub elems: Vec<Vec<usize>>,
    /// Boundary edges per patch: (node a, node b), a→b counterclockwise
    /// along the boundary.
    pub patches: Vec<(Patch, Vec<(usize, usize)>)>,
}

impl Mesh2 {
    /// Structured Q1 quads on `[0, lx] × [0, ly]`, `nx × ny` cells.
    #[must_use]
    pub fn quads(lx: f64, ly: f64, nx: usize, ny: usize) -> Mesh2 {
        Mesh2::mapped_quads(nx, ny, &|s, t| [lx * s, ly * t])
    }

    /// Structured P1 triangles (each grid cell split along its
    /// diagonal) on `[0, lx] × [0, ly]`.
    #[must_use]
    pub fn triangles(lx: f64, ly: f64, nx: usize, ny: usize) -> Mesh2 {
        let mut m = Mesh2::quads(lx, ly, nx, ny);
        let mut tris = Vec::with_capacity(2 * m.elems.len());
        for q in &m.elems {
            tris.push(vec![q[0], q[1], q[2]]);
            tris.push(vec![q[0], q[2], q[3]]);
        }
        m.elems = tris;
        m
    }

    /// Cook's membrane: the tapered panel with corners (0,0), (48,44),
    /// (48,60), (0,44), meshed with `n × n` mapped Q1 quads. Left edge
    /// clamps, right edge takes the shear load.
    #[must_use]
    pub fn cooks_membrane(n: usize) -> Mesh2 {
        Mesh2::mapped_quads(n, n, &|s, t| {
            let x = 48.0 * s;
            let y_bottom = 44.0 * s;
            let y_top = 44.0 + 16.0 * s;
            [x, y_bottom + (y_top - y_bottom) * t]
        })
    }

    /// Structured mapped Q1 quads over the unit parameter square.
    #[must_use]
    pub fn mapped_quads(nx: usize, ny: usize, map: &dyn Fn(f64, f64) -> [f64; 2]) -> Mesh2 {
        assert!(nx >= 1 && ny >= 1, "empty mesh");
        let id = |i: usize, j: usize| i + j * (nx + 1);
        let mut nodes = Vec::with_capacity((nx + 1) * (ny + 1));
        #[allow(clippy::cast_precision_loss)]
        for j in 0..=ny {
            for i in 0..=nx {
                nodes.push(map(i as f64 / nx as f64, j as f64 / ny as f64));
            }
        }
        let mut elems = Vec::with_capacity(nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                elems.push(vec![id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)]);
            }
        }
        let patches = vec![
            (
                Patch::Bottom,
                (0..nx).map(|i| (id(i, 0), id(i + 1, 0))).collect(),
            ),
            (
                Patch::Right,
                (0..ny).map(|j| (id(nx, j), id(nx, j + 1))).collect(),
            ),
            (
                Patch::Top,
                (0..nx).rev().map(|i| (id(i + 1, ny), id(i, ny))).collect(),
            ),
            (
                Patch::Left,
                (0..ny).rev().map(|j| (id(0, j + 1), id(0, j))).collect(),
            ),
        ];
        Mesh2 {
            nodes,
            elems,
            patches,
        }
    }

    /// Node count.
    #[must_use]
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// The edges of one patch.
    #[must_use]
    pub fn patch_edges(&self, p: Patch) -> Option<&[(usize, usize)]> {
        self.patches
            .iter()
            .find(|(q, _)| *q == p)
            .map(|(_, e)| e.as_slice())
    }

    /// All node ids lying on a patch.
    #[must_use]
    pub fn patch_nodes(&self, p: Patch) -> Vec<usize> {
        let mut ids: Vec<usize> = self
            .patch_edges(p)
            .map(|edges| edges.iter().flat_map(|&(a, b)| [a, b]).collect())
            .unwrap_or_default();
        ids.sort_unstable();
        ids.dedup();
        ids
    }
}

/// Shape values and reference-gradients for one element at a
/// quadrature point; returns (N, dN/dx, jacobian-weight).
#[must_use]
pub fn shapes_at(
    nodes: &[[f64; 2]],
    conn: &[usize],
    xi: f64,
    eta: f64,
) -> (Vec<f64>, Vec<[f64; 2]>, f64) {
    match conn.len() {
        3 => {
            // P1 triangle: reference (0,0)-(1,0)-(0,1).
            let n = vec![1.0 - xi - eta, xi, eta];
            let dref = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];
            iso_map(nodes, conn, &n, &dref)
        }
        4 => {
            // Q1 quad: reference [-1,1]².
            let n = vec![
                0.25 * (1.0 - xi) * (1.0 - eta),
                0.25 * (1.0 + xi) * (1.0 - eta),
                0.25 * (1.0 + xi) * (1.0 + eta),
                0.25 * (1.0 - xi) * (1.0 + eta),
            ];
            let dref = [
                [-0.25 * (1.0 - eta), -0.25 * (1.0 - xi)],
                [0.25 * (1.0 - eta), -0.25 * (1.0 + xi)],
                [0.25 * (1.0 + eta), 0.25 * (1.0 + xi)],
                [-0.25 * (1.0 + eta), 0.25 * (1.0 - xi)],
            ];
            iso_map(nodes, conn, &n, &dref)
        }
        k => unreachable!("element arity {k}"),
    }
}

fn iso_map(
    nodes: &[[f64; 2]],
    conn: &[usize],
    n: &[f64],
    dref: &[[f64; 2]],
) -> (Vec<f64>, Vec<[f64; 2]>, f64) {
    let mut j = [[0.0f64; 2]; 2];
    for (a, &node) in conn.iter().enumerate() {
        let p = nodes[node];
        for (r, jr) in j.iter_mut().enumerate() {
            jr[0] += dref[a][r] * p[0];
            jr[1] += dref[a][r] * p[1];
        }
    }
    let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
    let inv = [
        [j[1][1] / det, -j[0][1] / det],
        [-j[1][0] / det, j[0][0] / det],
    ];
    // Chain rule: dN/dx_c = sum_r dN/dxi_r * (J^-1)_{cr} — note the
    // TRANSPOSED inverse index (symmetric J on axis-aligned rectangles
    // hides a transpose slip; triangles and mapped panels do not).
    let grads = (0..conn.len())
        .map(|a| {
            [
                dref[a][0] * inv[0][0] + dref[a][1] * inv[0][1],
                dref[a][0] * inv[1][0] + dref[a][1] * inv[1][1],
            ]
        })
        .collect();
    (n.to_vec(), grads, det)
}

/// Quadrature points on the reference element: (xi, eta, weight).
#[must_use]
pub fn quad_points(arity: usize) -> Vec<(f64, f64, f64)> {
    if arity == 3 {
        // Degree-2 triangle rule (3 midpoints, weight 1/6 each on the
        // reference triangle of area 1/2).
        return vec![
            (0.5, 0.0, 1.0 / 6.0),
            (0.5, 0.5, 1.0 / 6.0),
            (0.0, 0.5, 1.0 / 6.0),
        ];
    }
    // 2×2 Gauss on [-1,1]².
    let g = 1.0 / 3.0f64.sqrt();
    vec![(-g, -g, 1.0), (g, -g, 1.0), (g, g, 1.0), (-g, g, 1.0)]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quad_mesh_area_is_exact() {
        let m = Mesh2::quads(2.0, 1.0, 4, 3);
        let mut area = 0.0;
        for e in &m.elems {
            for &(xi, eta, w) in &quad_points(4) {
                let (_, _, det) = shapes_at(&m.nodes, e, xi, eta);
                area += w * det;
            }
        }
        assert!((area - 2.0).abs() < 1e-12);
    }

    #[test]
    fn cooks_area_matches_trapezoid() {
        let m = Mesh2::cooks_membrane(8);
        let mut area = 0.0;
        for e in &m.elems {
            for &(xi, eta, w) in &quad_points(4) {
                let (_, _, det) = shapes_at(&m.nodes, e, xi, eta);
                area += w * det;
            }
        }
        // Trapezoid: 48 × (44 + 16)/2.
        assert!((area - 48.0 * 30.0).abs() < 1e-9, "area {area}");
    }
}

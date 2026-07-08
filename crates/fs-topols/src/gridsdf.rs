//! The discrete level set: nodal φ on a uniform grid, bilinear between
//! nodes, implementing fs-cutfem's [`CutSdf`] so THE LEVEL SET IS THE
//! GEOMETRY — physics evaluates directly on the evolving field with
//! zero meshing.
//!
//! The enclosure is EXACT (up to an outward roundoff pad): a bilinear
//! restricted to any axis-aligned rectangle attains its extrema at the
//! rectangle's corners (it is affine along every axis-aligned line),
//! so a box's range is the hull of per-grid-cell clipped-corner
//! evaluations — the certified-classification contract fs-cutfem
//! requires.

use fs_cutfem::CutSdf;
use fs_ivl::Interval;

/// A nodal level-set field on the `(n+1)²` lattice over `[0,1]²`.
#[derive(Debug, Clone)]
pub struct GridSdf {
    n: usize,
    phi: Vec<f64>,
}

impl GridSdf {
    /// Sample a function at the nodes.
    #[must_use]
    pub fn from_fn(n: usize, f: &dyn Fn(f64, f64) -> f64) -> GridSdf {
        assert!(n >= 2, "grid too small");
        let mut phi = Vec::with_capacity((n + 1) * (n + 1));
        #[allow(clippy::cast_precision_loss)]
        for j in 0..=n {
            for i in 0..=n {
                phi.push(f(i as f64 / n as f64, j as f64 / n as f64));
            }
        }
        GridSdf { n, phi }
    }

    /// Cells per side.
    #[must_use]
    pub fn n(&self) -> usize {
        self.n
    }

    /// Grid spacing.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn h(&self) -> f64 {
        1.0 / self.n as f64
    }

    /// Nodal value.
    #[must_use]
    pub fn node(&self, i: usize, j: usize) -> f64 {
        self.phi[i + j * (self.n + 1)]
    }

    /// Mutable nodal value.
    pub fn node_mut(&mut self, i: usize, j: usize) -> &mut f64 {
        &mut self.phi[i + j * (self.n + 1)]
    }

    /// The raw nodal slice (row-major, x fastest).
    #[must_use]
    pub fn nodes(&self) -> &[f64] {
        &self.phi
    }

    /// Mutable raw nodal slice.
    pub fn nodes_mut(&mut self) -> &mut [f64] {
        &mut self.phi
    }

    /// Node position.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn pos(&self, i: usize, j: usize) -> [f64; 2] {
        [i as f64 / self.n as f64, j as f64 / self.n as f64]
    }

    /// Clamp a point into the grid and find its cell + local coords.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    fn locate(&self, p: [f64; 2]) -> (usize, usize, f64, f64) {
        #[allow(clippy::cast_precision_loss)]
        let nf = self.n as f64;
        let x = (p[0] * nf).clamp(0.0, nf - 1e-12);
        let y = (p[1] * nf).clamp(0.0, nf - 1e-12);
        let ci = (x.floor() as usize).min(self.n - 1);
        let cj = (y.floor() as usize).min(self.n - 1);
        #[allow(clippy::cast_precision_loss)]
        let (xi, et) = (x - ci as f64, y - cj as f64);
        (ci, cj, xi, et)
    }

    /// Bilinear value at a point.
    #[must_use]
    pub fn value_at(&self, p: [f64; 2]) -> f64 {
        let (ci, cj, xi, et) = self.locate(p);
        let v00 = self.node(ci, cj);
        let v10 = self.node(ci + 1, cj);
        let v01 = self.node(ci, cj + 1);
        let v11 = self.node(ci + 1, cj + 1);
        v00 * (1.0 - xi) * (1.0 - et)
            + v10 * xi * (1.0 - et)
            + v01 * (1.0 - xi) * et
            + v11 * xi * et
    }

    /// Bilinear gradient at a point.
    #[must_use]
    pub fn gradient_at(&self, p: [f64; 2]) -> [f64; 2] {
        let (ci, cj, xi, et) = self.locate(p);
        let v00 = self.node(ci, cj);
        let v10 = self.node(ci + 1, cj);
        let v01 = self.node(ci, cj + 1);
        let v11 = self.node(ci + 1, cj + 1);
        let inv_h = 1.0 / self.h();
        [
            ((v10 - v00) * (1.0 - et) + (v11 - v01) * et) * inv_h,
            ((v01 - v00) * (1.0 - xi) + (v11 - v10) * xi) * inv_h,
        ]
    }
}

impl CutSdf for GridSdf {
    fn value(&self, p: [f64; 2]) -> f64 {
        self.value_at(p)
    }

    fn gradient(&self, p: [f64; 2]) -> [f64; 2] {
        self.gradient_at(p)
    }

    fn enclose(&self, lo: [f64; 2], hi: [f64; 2]) -> Interval {
        // Hull of clipped-corner evaluations over every overlapped
        // grid cell (bilinear extrema sit at rectangle corners).
        #[allow(clippy::cast_precision_loss)]
        let nf = self.n as f64;
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let ci0 = ((lo[0] * nf).floor().clamp(0.0, nf - 1.0)) as usize;
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let ci1 = ((hi[0] * nf).ceil().clamp(1.0, nf)) as usize;
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let cj0 = ((lo[1] * nf).floor().clamp(0.0, nf - 1.0)) as usize;
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let cj1 = ((hi[1] * nf).ceil().clamp(1.0, nf)) as usize;
        let mut vmin = f64::INFINITY;
        let mut vmax = f64::NEG_INFINITY;
        let h = self.h();
        for cj in cj0..cj1 {
            for ci in ci0..ci1 {
                #[allow(clippy::cast_precision_loss)]
                let (cx0, cy0) = (ci as f64 * h, cj as f64 * h);
                let (cx1, cy1) = (cx0 + h, cy0 + h);
                let xa = lo[0].max(cx0);
                let xb = hi[0].min(cx1);
                let ya = lo[1].max(cy0);
                let yb = hi[1].min(cy1);
                if xa > xb || ya > yb {
                    continue;
                }
                for p in [[xa, ya], [xb, ya], [xb, yb], [xa, yb]] {
                    let v = self.value_at(p);
                    vmin = vmin.min(v);
                    vmax = vmax.max(v);
                }
            }
        }
        if !vmin.is_finite() || !vmax.is_finite() {
            return Interval::WHOLE;
        }
        // Outward roundoff pad (a few ulps of the magnitude scale).
        let pad = 1e-13 * (1.0 + vmin.abs().max(vmax.abs()));
        Interval::new(vmin - pad, vmax + pad)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn enclosure_contains_dense_samples() {
        let g = GridSdf::from_fn(16, &|x, y| {
            ((x - 0.5) * (x - 0.5) + (y - 0.4) * (y - 0.4)).sqrt() - 0.3
        });
        let boxes = [
            ([0.1, 0.1], [0.3, 0.35]),
            ([0.45, 0.45], [0.55, 0.55]),
            ([0.0, 0.7], [1.0, 1.0]),
        ];
        for (lo, hi) in boxes {
            let iv = g.enclose(lo, hi);
            for si in 0..=20 {
                for sj in 0..=20 {
                    let p = [
                        lo[0] + (hi[0] - lo[0]) * f64::from(si) / 20.0,
                        lo[1] + (hi[1] - lo[1]) * f64::from(sj) / 20.0,
                    ];
                    assert!(iv.contains(g.value_at(p)), "containment law");
                }
            }
        }
    }
}

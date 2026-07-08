//! DE-HOMOGENIZATION (bead 7tv.14): the graded density field realized
//! as EXPLICIT micro-geometry — one hole per macro cell, radius from
//! inverting the hole-area/density relation r = h·√((1−ρ)/π) — and
//! re-analyzed at FULL resolution with fs-cutfem CutElasticity on the
//! hole-array SDF (traction-free hole boundaries, ghost-stabilized cut
//! cells, no meshing). The verification gates live in the battery:
//! the full-resolution compliance sits inside an HONESTY BAND of the
//! homogenized prediction, and — the design-transfer statement that
//! actually matters — the graded hole array still beats the uniform
//! hole array after realization. Load conventions match the macro
//! model exactly (left clamp, uniform downward right-edge traction,
//! trapezoidal edge quadrature).

use fs_cutfem::{CutSdf, Quadtree};
use fs_solid::SolidError;
use fs_solid::cutfront::CutElasticity;

/// The realized micro-geometry: one circular hole per macro cell.
pub struct HoleArray {
    /// Hole centers.
    pub centers: Vec<[f64; 2]>,
    /// Hole radii.
    pub radii: Vec<f64>,
}

impl HoleArray {
    /// Realize a graded density field on an `nx × ny` unit-square
    /// macro grid: hole area per cell = (1 − ρ)·cell area. Holes
    /// below `min_r` are dropped (a near-solid cell is solid — the
    /// cut solver cannot resolve point holes, and physically they do
    /// not matter).
    #[must_use]
    pub fn realize(rho: &[f64], nx: usize, ny: usize, min_r: f64) -> HoleArray {
        let (hx, hy) = (1.0 / nx as f64, 1.0 / ny as f64);
        let mut centers = Vec::new();
        let mut radii = Vec::new();
        for j in 0..ny {
            for i in 0..nx {
                let dens = rho[j * nx + i].clamp(0.0, 1.0);
                let r = (hx * hy * (1.0 - dens) / core::f64::consts::PI).sqrt();
                if r >= min_r {
                    centers.push([(i as f64 + 0.5) * hx, (j as f64 + 0.5) * hy]);
                    radii.push(r);
                }
            }
        }
        HoleArray { centers, radii }
    }

    /// Smallest and largest realized radius (0, 0) when empty.
    #[must_use]
    pub fn radius_range(&self) -> (f64, f64) {
        if self.radii.is_empty() {
            return (0.0, 0.0);
        }
        self.radii
            .iter()
            .fold((f64::INFINITY, 0.0f64), |(lo, hi), &r| {
                (lo.min(r), hi.max(r))
            })
    }
}

impl CutSdf for HoleArray {
    /// Material = complement of the holes: φ = max_k (r_k − |x − c_k|)
    /// (negative in material, positive inside a hole).
    fn value(&self, p: [f64; 2]) -> f64 {
        if self.radii.is_empty() {
            return -1.0;
        }
        let mut phi = f64::NEG_INFINITY;
        for (c, &r) in self.centers.iter().zip(&self.radii) {
            let d = (p[0] - c[0]).hypot(p[1] - c[1]);
            phi = phi.max(r - d);
        }
        phi
    }

    fn gradient(&self, p: [f64; 2]) -> [f64; 2] {
        // Gradient of the winning hole's term: −(x − c)/|x − c|.
        let mut best = f64::NEG_INFINITY;
        let mut g = [0.0f64, 0.0];
        for (c, &r) in self.centers.iter().zip(&self.radii) {
            let dx = p[0] - c[0];
            let dy = p[1] - c[1];
            let d = dx.hypot(dy).max(1e-30);
            let phi = r - d;
            if phi > best {
                best = phi;
                g = [-dx / d, -dy / d];
            }
        }
        g
    }

    fn enclose(&self, lo: [f64; 2], hi: [f64; 2]) -> fs_ivl::Interval {
        if self.radii.is_empty() {
            return fs_ivl::Interval::point(-1.0);
        }
        // Per hole: |x − c| over the box lies in [dmin, dmax] with
        // dmin the box-point distance and dmax the farthest corner;
        // φ_k ∈ [r − dmax, r − dmin]. max over holes of enclosures
        // encloses the max (monotone).
        let mut out_lo = f64::NEG_INFINITY;
        let mut out_hi = f64::NEG_INFINITY;
        for (c, &r) in self.centers.iter().zip(&self.radii) {
            let cx = c[0].clamp(lo[0], hi[0]);
            let cy = c[1].clamp(lo[1], hi[1]);
            let dmin = (cx - c[0]).hypot(cy - c[1]);
            let mut dmax = 0.0f64;
            for corner in [
                [lo[0], lo[1]],
                [hi[0], lo[1]],
                [lo[0], hi[1]],
                [hi[0], hi[1]],
            ] {
                dmax = dmax.max((corner[0] - c[0]).hypot(corner[1] - c[1]));
            }
            out_lo = out_lo.max(r - dmax);
            out_hi = out_hi.max(r - dmin);
        }
        fs_ivl::Interval::new(out_lo, out_hi)
    }
}

/// Full-resolution CutFEM compliance of a hole array under the macro
/// load convention: left edge clamped, uniform downward traction of
/// magnitude 1 on the right edge, compliance by trapezoidal edge
/// quadrature of t·u.
///
/// # Errors
/// Propagates [`SolidError`] from the cut solver.
pub fn fullres_compliance(holes: &HoleArray, level: u32) -> Result<f64, SolidError> {
    let grid = Quadtree::uniform(level);
    let clamp = |x: f64, _y: f64| x < 1e-9;
    let traction = |x: f64, _y: f64| -> [f64; 2] {
        if x > 1.0 - 1e-9 {
            [0.0, -1.0]
        } else {
            [0.0, 0.0]
        }
    };
    let solver = CutElasticity {
        grid: &grid,
        sdf: holes,
        youngs: 2.6,
        poisson: 0.3,
        nitsche_beta: 20.0,
        ghost_gamma: 0.5,
        quad_depth: 2,
        clamp: Some(&clamp),
        boundary_traction: Some(&traction),
        traction_free_interface: true,
    };
    let sol = solver.solve(&|_, _| [0.0, 0.0], &|_, _| [0.0, 0.0])?;
    let n = 1usize << level;
    let h = 1.0 / n as f64;
    let mut compliance = 0.0;
    for (&(gi, gj), u) in sol.nodal() {
        if gi as usize == n {
            let w = if gj == 0 || gj as usize == n {
                0.5 * h
            } else {
                h
            };
            compliance -= w * u[1];
        }
    }
    Ok(compliance)
}

//! fs-viz — scientific visualization primitives. Layer: L5.
//!
//! A 10⁸-cell field is illegible as raw numbers; the job of visualization is to
//! compress it into a picture an agent can read in one glance. This v0 is the
//! verifiable core of that job — the pieces with ANALYTIC ground truth, so a
//! rendered topology is a certificate rather than a guess:
//!
//! - [`streamline`] — RK4 integration of a flow field; on a rotation field the
//!   streamline is a circle (radius conserved), on a saddle field it is a
//!   hyperbola (`xy` conserved);
//! - [`classify_hessian`] — the Morse critical-point type (min / max / saddle)
//!   and index from a scalar field's Hessian — the atom of a Morse–Smale
//!   complex;
//! - [`Grid2::isocontour_crossings`] — the isocontour edge crossings of a scalar
//!   grid (the one contouring implementation), which on a circle SDF all lie on
//!   the circle.
//! - [`Grid3::isosurface`] — bounded deterministic marching tetrahedra from an
//!   owned x-fastest scalar grid to a renderer-ready indexed triangle mesh.
//! - [`ScalarField3`] — a bounded versioned artifact codec with explicit
//!   node/cell centering, quantity, and units for ledger composition.
//!
//! Deterministic; no dependencies.

mod isosurface;
mod scalar_field;

pub use isosurface::{Grid3, Grid3Error, IsoMesh3, IsoSurfaceError, Vec3};
pub use scalar_field::{
    SCALAR_FIELD3_ARTIFACT_KIND, SCALAR_FIELD3_SCHEMA_VERSION, ScalarField3, ScalarField3Error,
    ScalarFieldSemantics, ScalarLayout3,
};

/// A 2-D point / vector.
pub type Vec2 = [f64; 2];

/// Integrate a streamline of `field` from `seed` with `steps` RK4 steps of size
/// `dt`. Returns the ordered polyline (including the seed).
#[must_use]
pub fn streamline(field: impl Fn(Vec2) -> Vec2, seed: Vec2, dt: f64, steps: usize) -> Vec<Vec2> {
    let mut pts = Vec::with_capacity(steps + 1);
    pts.push(seed);
    let mut p = seed;
    for _ in 0..steps {
        let k1 = field(p);
        let k2 = field([p[0] + 0.5 * dt * k1[0], p[1] + 0.5 * dt * k1[1]]);
        let k3 = field([p[0] + 0.5 * dt * k2[0], p[1] + 0.5 * dt * k2[1]]);
        let k4 = field([p[0] + dt * k3[0], p[1] + dt * k3[1]]);
        p = [
            p[0] + dt / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]),
            p[1] + dt / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]),
        ];
        pts.push(p);
    }
    pts
}

/// The Morse type of a critical point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CriticalKind {
    /// A local minimum (Morse index 0).
    Minimum,
    /// A saddle (Morse index 1 in 2-D).
    Saddle,
    /// A local maximum (Morse index 2 in 2-D).
    Maximum,
    /// Degenerate (a zero Hessian eigenvalue) — not a Morse critical point.
    Degenerate,
}

/// A classified critical point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CriticalPoint {
    /// The Morse type.
    pub kind: CriticalKind,
    /// The Morse index (number of negative Hessian eigenvalues).
    pub morse_index: usize,
}

/// Classify a critical point from its (symmetric) `2×2` Hessian: the Morse index
/// is the number of negative eigenvalues; a near-zero eigenvalue is degenerate.
#[must_use]
pub fn classify_hessian(hessian: [[f64; 2]; 2], tol: f64) -> CriticalPoint {
    let (a, b, c) = (hessian[0][0], hessian[0][1], hessian[1][1]);
    let mean = f64::midpoint(a, c);
    let disc = (((a - c) / 2.0).powi(2) + b * b).sqrt();
    let (l1, l2) = (mean - disc, mean + disc);
    if l1.abs() <= tol || l2.abs() <= tol {
        return CriticalPoint {
            kind: CriticalKind::Degenerate,
            morse_index: usize::from(l1 < -tol) + usize::from(l2 < -tol),
        };
    }
    let morse_index = usize::from(l1 < 0.0) + usize::from(l2 < 0.0);
    let kind = match morse_index {
        0 => CriticalKind::Minimum,
        1 => CriticalKind::Saddle,
        _ => CriticalKind::Maximum,
    };
    CriticalPoint { kind, morse_index }
}

/// A scalar field sampled on a regular 2-D grid.
#[derive(Debug, Clone)]
pub struct Grid2 {
    nx: usize,
    ny: usize,
    lo: Vec2,
    hi: Vec2,
    values: Vec<f64>,
}

impl Grid2 {
    /// Sample `f` on an `nx × ny` regular grid spanning `[lo, hi]`.
    ///
    /// # Panics
    /// If `nx < 2` or `ny < 2`.
    #[must_use]
    pub fn from_fn(nx: usize, ny: usize, lo: Vec2, hi: Vec2, f: impl Fn(Vec2) -> f64) -> Grid2 {
        assert!(nx >= 2 && ny >= 2, "grid needs at least 2 points per axis");
        let mut values = Vec::with_capacity(nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                let p = grid_point(nx, ny, lo, hi, i, j);
                values.push(f(p));
            }
        }
        Grid2 {
            nx,
            ny,
            lo,
            hi,
            values,
        }
    }

    /// The scalar value at grid node `(i, j)`.
    #[must_use]
    pub fn at(&self, i: usize, j: usize) -> f64 {
        self.values[j * self.nx + i]
    }

    /// The world coordinate of grid node `(i, j)`.
    #[must_use]
    pub fn point(&self, i: usize, j: usize) -> Vec2 {
        grid_point(self.nx, self.ny, self.lo, self.hi, i, j)
    }

    /// The isocontour edge crossings at `iso`: for every grid edge whose two
    /// endpoints straddle `iso`, the linearly-interpolated crossing point (the
    /// one contouring pass shared with SDF→mesh conversion).
    #[must_use]
    pub fn isocontour_crossings(&self, iso: f64) -> Vec<Vec2> {
        let mut out = Vec::new();
        let mut cross = |a: Vec2, va: f64, b: Vec2, vb: f64| {
            let (da, db) = (va - iso, vb - iso);
            if (da < 0.0) != (db < 0.0) && (da != 0.0 || db != 0.0) {
                let t = da / (da - db); // da + t(db-da) = 0
                out.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]);
            }
        };
        for j in 0..self.ny {
            for i in 0..self.nx {
                let (p, v) = (self.point(i, j), self.at(i, j));
                if i + 1 < self.nx {
                    cross(p, v, self.point(i + 1, j), self.at(i + 1, j));
                }
                if j + 1 < self.ny {
                    cross(p, v, self.point(i, j + 1), self.at(i, j + 1));
                }
            }
        }
        out
    }
}

fn grid_point(nx: usize, ny: usize, lo: Vec2, hi: Vec2, i: usize, j: usize) -> Vec2 {
    [
        lo[0] + (hi[0] - lo[0]) * i as f64 / (nx - 1) as f64,
        lo[1] + (hi[1] - lo[1]) * j as f64 / (ny - 1) as f64,
    ]
}

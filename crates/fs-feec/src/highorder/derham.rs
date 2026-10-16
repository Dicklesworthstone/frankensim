//! The full tensor-product de Rham complex on structured hex grids
//! (tfz.6 slice 2): H¹ → H(curl) → H(div) → L² with EXACT discrete
//! derivatives.
//!
//! 1D structure: the continuous order-r Lobatto space C_r (slice 1's
//! lattice) and the discontinuous order-(r−1) LEGENDRE space D_{r−1}
//! (r modes per cell, orthogonal ⇒ diagonal 1D mass). The derivative
//! of the Lobatto basis lands EXACTLY in the Legendre basis
//! (N'_k = √((2k−1)/2)·L_{k−1} — the integrated-Legendre identity),
//! so the 1D derivative operator G: C_r → D_{r−1} has closed-form
//! entries and the 3D operators are Kronecker assemblies:
//!
//! - grad: (G⊗I⊗I, I⊗G⊗I, I⊗I⊗G)
//! - curl_x = (I⊗G⊗I)·u_z − (I⊗I⊗G)·u_y (cyclic)
//! - div = Σ ∂_axis on the matching component
//!
//! curl∘grad and div∘curl vanish to ROUNDOFF: the two paths of each
//! mixed second derivative multiply the same two G entries in
//! different association orders, so cancellation is at the ε level
//! (the exact-integer dd = 0 of the lowest-order complex trades into
//! floating-point commutation here; the battery pins the bound).
//!
//! Component factor types: E = ((D,C,C),(C,D,C),(C,C,D)),
//! F = ((C,D,D),(D,C,D),(D,D,C)), W = (D,D,D).

use crate::highorder::hex::TensorSpace;
use crate::highorder::quad1d::{gauss_legendre, legendre, lobatto_shapes};
use fs_math::det;

/// A dense rectangular 1D operator (rows × cols, row-major).
pub struct Mat1 {
    /// Row count.
    pub rows: usize,
    /// Column count.
    pub cols: usize,
    /// Row-major entries.
    pub a: Vec<f64>,
}

/// The tensor-product de Rham complex on an m³ grid of the unit cube
/// at order r (C_r / D_{r−1} factor pair).
pub struct TensorDeRham {
    /// Cells per axis.
    pub m: usize,
    /// Order r ≥ 1.
    pub r: usize,
    /// Continuous 1D dof count m·r + 1.
    pub nc: usize,
    /// Discontinuous 1D dof count m·r.
    pub nd: usize,
    /// 1D derivative operator G (nd × nc), closed form.
    pub g: Mat1,
    /// Assembled continuous 1D mass (nc × nc dense).
    pub mass_c: Vec<f64>,
    /// Diagonal of the discontinuous 1D mass (Legendre orthogonality).
    pub mass_d: Vec<f64>,
}

impl TensorDeRham {
    /// Build the complex.
    ///
    /// # Panics
    /// If `m == 0` or `r == 0`.
    #[must_use]
    pub fn new(m: usize, r: usize) -> TensorDeRham {
        assert!(m >= 1 && r >= 1, "TensorDeRham needs m >= 1, r >= 1");
        let h = 1.0 / m as f64;
        let sp = TensorSpace::new(m, r);
        let (mass_c, _) = sp.assembled_1d();
        let (nc, nd) = (m * r + 1, m * r);
        // G: cell c, Legendre row c·r + j. Physical d/dx = (2/h)·d/dξ.
        // Vertex chain: N'_0 = −1/h, N'_1 = +1/h into L_0 of the cell;
        // bubble k: N'_k = (2/h)·√((2k−1)/2)·L_{k−1}.
        let mut g = vec![0.0f64; nd * nc];
        for c in 0..m {
            let left = c * r;
            let right = (c + 1) * r;
            g[(c * r) * nc + left] = -1.0 / h;
            g[(c * r) * nc + right] = 1.0 / h;
            for k in 2..=r {
                let row = c * r + (k - 1);
                let col = c * r + (k - 1);
                g[row * nc + col] = (2.0 / h) * det::sqrt((2.0 * k as f64 - 1.0) / 2.0);
            }
        }
        // Legendre 1D mass diagonal: ∫_cell L_j² = (h/2)·2/(2j+1).
        let mut mass_d = vec![0.0f64; nd];
        for c in 0..m {
            for j in 0..r {
                mass_d[c * r + j] = (h / 2.0) * 2.0 / (2.0 * j as f64 + 1.0);
            }
        }
        TensorDeRham {
            m,
            r,
            nc,
            nd,
            g: Mat1 {
                rows: nd,
                cols: nc,
                a: g,
            },
            mass_c,
            mass_d,
        }
    }

    /// Dof counts (S, E, F, W).
    #[must_use]
    pub fn dims(&self) -> [usize; 4] {
        let (c, d) = (self.nc, self.nd);
        [c * c * c, 3 * d * c * c, 3 * c * d * d, d * d * d]
    }

    /// grad: S → E (component blocks x, y, z concatenated).
    #[must_use]
    pub fn grad(&self, s: &[f64]) -> Vec<f64> {
        let (c, d) = (self.nc, self.nd);
        assert_eq!(s.len(), c * c * c, "grad input length");
        let mut out = Vec::with_capacity(3 * d * c * c);
        out.extend(apply_axis(&self.g, s, [c, c, c], 0));
        out.extend(apply_axis(&self.g, s, [c, c, c], 1));
        out.extend(apply_axis(&self.g, s, [c, c, c], 2));
        out
    }

    /// curl: E → F.
    #[must_use]
    pub fn curl(&self, e: &[f64]) -> Vec<f64> {
        let (c, d) = (self.nc, self.nd);
        let bx = d * c * c;
        assert_eq!(e.len(), 3 * bx, "curl input length");
        let (ex, rest) = e.split_at(bx);
        let (ey, ez) = rest.split_at(bx);
        // F_x = ∂_y e_z − ∂_z e_y : e_z ∈ (C,C,D), e_y ∈ (C,D,C).
        let mut fx = apply_axis(&self.g, ez, [c, c, d], 1);
        for (a, b) in fx.iter_mut().zip(apply_axis(&self.g, ey, [c, d, c], 2)) {
            *a -= b;
        }
        // F_y = ∂_z e_x − ∂_x e_z : e_x ∈ (D,C,C), e_z ∈ (C,C,D).
        let mut fy = apply_axis(&self.g, ex, [d, c, c], 2);
        for (a, b) in fy.iter_mut().zip(apply_axis(&self.g, ez, [c, c, d], 0)) {
            *a -= b;
        }
        // F_z = ∂_x e_y − ∂_y e_x : e_y ∈ (C,D,C), e_x ∈ (D,C,C).
        let mut fz = apply_axis(&self.g, ey, [c, d, c], 0);
        for (a, b) in fz.iter_mut().zip(apply_axis(&self.g, ex, [d, c, c], 1)) {
            *a -= b;
        }
        let mut out = fx;
        out.extend(fy);
        out.extend(fz);
        out
    }

    /// div: F → W.
    #[must_use]
    pub fn div(&self, f: &[f64]) -> Vec<f64> {
        let (c, d) = (self.nc, self.nd);
        let bx = c * d * d;
        assert_eq!(f.len(), 3 * bx, "div input length");
        let (fx, rest) = f.split_at(bx);
        let (fy, fz) = rest.split_at(bx);
        let mut out = apply_axis(&self.g, fx, [c, d, d], 0);
        for (a, b) in out.iter_mut().zip(apply_axis(&self.g, fy, [d, c, d], 1)) {
            *a += b;
        }
        for (a, b) in out.iter_mut().zip(apply_axis(&self.g, fz, [d, d, c], 2)) {
            *a += b;
        }
        out
    }

    /// 1D canonical projection dofs: continuous π_C of a scalar
    /// function (endpoint values + bubble dofs from Legendre moments
    /// of f′ — the construction that makes d∘π_C = π_D∘d hold by
    /// design), evaluated per axis on the [0,1] mesh.
    #[must_use]
    pub fn project_c_1d<F: Fn(f64) -> f64, D: Fn(f64) -> f64>(&self, f: &F, df: &D) -> Vec<f64> {
        let h = 1.0 / self.m as f64;
        let mut dofs = vec![0.0f64; self.nc];
        let (qx, qw) = gauss_legendre(self.r + 4);
        for c in 0..self.m {
            let (xl, xr) = (c as f64 * h, (c + 1) as f64 * h);
            dofs[c * self.r] = f(xl);
            dofs[(c + 1) * self.r] = f(xr);
            // Bubble k dof = c_{k−1}/√((2k−1)/2) with c_j the Legendre
            // coefficient of the REFERENCE derivative (h/2)·f′.
            for k in 2..=self.r {
                let j = k - 1;
                let mut moment = 0.0f64;
                for (&xi, &wi) in qx.iter().zip(&qw) {
                    let x = f64::midpoint(xl, xr) + xi * h / 2.0;
                    let (lj, _) = legendre(j, xi);
                    moment += wi * df(x) * (h / 2.0) * lj;
                }
                let cj = (j as f64 + 0.5) * moment;
                dofs[c * self.r + k - 1] = cj / det::sqrt(j as f64 + 0.5);
            }
        }
        dofs
    }

    /// 1D discontinuous π_D dofs: per-cell Legendre coefficients of f.
    #[must_use]
    pub fn project_d_1d<F: Fn(f64) -> f64>(&self, f: &F) -> Vec<f64> {
        let h = 1.0 / self.m as f64;
        let mut dofs = vec![0.0f64; self.nd];
        let (qx, qw) = gauss_legendre(self.r + 4);
        for c in 0..self.m {
            let xm = (c as f64).mul_add(h, h / 2.0);
            for j in 0..self.r {
                let mut moment = 0.0f64;
                for (&xi, &wi) in qx.iter().zip(&qw) {
                    let (lj, _) = legendre(j, xi);
                    moment += wi * f(xm + xi * h / 2.0) * lj;
                }
                dofs[c * self.r + j] = (j as f64 + 0.5) * moment;
            }
        }
        dofs
    }

    /// Evaluate a 1D continuous dof vector at physical point x.
    #[must_use]
    pub fn eval_c_1d(&self, dofs: &[f64], x: f64) -> f64 {
        let h = 1.0 / self.m as f64;
        let c = ((x / h) as usize).min(self.m - 1);
        let xi = ((x - c as f64 * h) / h).mul_add(2.0, -1.0);
        let (vals, _) = lobatto_shapes(self.r, xi);
        let mut acc = 0.0f64;
        for (l, v) in vals.iter().enumerate() {
            let g = match l {
                0 => c * self.r,
                1 => (c + 1) * self.r,
                k => c * self.r + k - 1,
            };
            acc = v.mul_add(dofs[g], acc);
        }
        acc
    }

    /// Evaluate a 1D discontinuous dof vector at physical point x.
    #[must_use]
    pub fn eval_d_1d(&self, dofs: &[f64], x: f64) -> f64 {
        let h = 1.0 / self.m as f64;
        let c = ((x / h) as usize).min(self.m - 1);
        let xi = ((x - c as f64 * h) / h).mul_add(2.0, -1.0);
        let mut acc = 0.0f64;
        for j in 0..self.r {
            let (lj, _) = legendre(j, xi);
            acc = lj.mul_add(dofs[c * self.r + j], acc);
        }
        acc
    }

    /// L2 error of a 1D discontinuous dof vector against an analytic f.
    #[must_use]
    pub fn l2_error_d_1d<F: Fn(f64) -> f64>(&self, dofs: &[f64], f: &F) -> f64 {
        let h = 1.0 / self.m as f64;
        let (qx, qw) = gauss_legendre(self.r + 4);
        let mut total = 0.0f64;
        for c in 0..self.m {
            let xm = (c as f64).mul_add(h, h / 2.0);
            for (&xi, &wi) in qx.iter().zip(&qw) {
                let x = xm + xi * h / 2.0;
                let e = self.eval_d_1d(dofs, x) - f(x);
                total += wi * (h / 2.0) * e * e;
            }
        }
        det::sqrt(total)
    }

    /// L2 error of a 1D continuous dof vector against an analytic f.
    #[must_use]
    pub fn l2_error_c_1d<F: Fn(f64) -> f64>(&self, dofs: &[f64], f: &F) -> f64 {
        let h = 1.0 / self.m as f64;
        let (qx, qw) = gauss_legendre(self.r + 4);
        let mut total = 0.0f64;
        for c in 0..self.m {
            let xm = (c as f64).mul_add(h, h / 2.0);
            for (&xi, &wi) in qx.iter().zip(&qw) {
                let x = xm + xi * h / 2.0;
                let e = self.eval_c_1d(dofs, x) - f(x);
                total += wi * (h / 2.0) * e * e;
            }
        }
        det::sqrt(total)
    }
}

/// Apply a rectangular 1D operator along one axis of a flat 3D array
/// with per-axis extents `shape` (row-major x-major layout, matching
/// `TensorSpace::gid`). Returns the transformed array (the chosen
/// axis's extent becomes `op.rows`).
#[must_use]
pub fn apply_axis(op: &Mat1, src: &[f64], shape: [usize; 3], axis: usize) -> Vec<f64> {
    let [sx, sy, sz] = shape;
    assert_eq!(src.len(), sx * sy * sz, "apply_axis shape mismatch");
    assert_eq!(op.cols, shape[axis], "apply_axis operator/axis mismatch");
    let mut out_shape = shape;
    out_shape[axis] = op.rows;
    let [ox, oy, oz] = out_shape;
    let mut out = vec![0.0f64; ox * oy * oz];
    match axis {
        0 => {
            for i in 0..ox {
                for l in 0..sx {
                    let ail = op.a[i * op.cols + l];
                    if ail != 0.0 {
                        for j in 0..sy {
                            for k in 0..sz {
                                out[(i * oy + j) * oz + k] = ail.mul_add(
                                    src[(l * sy + j) * sz + k],
                                    out[(i * oy + j) * oz + k],
                                );
                            }
                        }
                    }
                }
            }
        }
        1 => {
            for i in 0..ox {
                for j in 0..oy {
                    for l in 0..sy {
                        let ajl = op.a[j * op.cols + l];
                        if ajl != 0.0 {
                            for k in 0..sz {
                                out[(i * oy + j) * oz + k] = ajl.mul_add(
                                    src[(i * sy + l) * sz + k],
                                    out[(i * oy + j) * oz + k],
                                );
                            }
                        }
                    }
                }
            }
        }
        2 => {
            for i in 0..ox {
                for j in 0..oy {
                    for k in 0..oz {
                        let mut acc = 0.0f64;
                        for l in 0..sz {
                            acc = op.a[k * op.cols + l].mul_add(src[(i * sy + j) * sz + l], acc);
                        }
                        out[(i * oy + j) * oz + k] = acc;
                    }
                }
            }
        }
        _ => panic!("axis must be 0..=2"),
    }
    out
}

//! Tensor-product H¹ space Q_r on structured hex grids of the unit
//! cube, with SUM-FACTORIZED matrix-free apply (tfz.6 slice 1).
//!
//! The 1D global dof lattice per axis has n₁ = m·r + 1 dofs: cell c
//! owns [c·r, c·r + r], vertex-chain dofs sit at multiples of r, and
//! Lobatto bubbles fill between (a bubble vanishes at both cell ends,
//! so only lattice dofs 0 and m·r have nonzero trace on the domain
//! faces — the Dirichlet logic). The 3D space is the tensor of three
//! lattices; the global Galerkin Poisson operator IS
//! K₁⊗M₁⊗M₁ + M₁⊗K₁⊗M₁ + M₁⊗M₁⊗K₁ over ASSEMBLED 1D operators
//! (basis separability + assembly linearity), which the battery uses
//! as the roundoff-match reference.
//!
//! The apply path is the P6 doctrine: per element, gather (r+1)³
//! local dofs, contract dimension by dimension with the dense
//! (r+1)² 1D element matrices — O(r⁴) per element instead of the
//! naive O(r⁶) — and scatter-add in fixed element order (bitwise
//! deterministic run-to-run and across ISAs).

use crate::highorder::quad1d::{element_matrices, gauss_legendre, lobatto_shapes};

/// A Q_r tensor-product H¹ space on an m×m×m grid of the unit cube.
pub struct TensorSpace {
    /// Cells per axis.
    pub m: usize,
    /// Polynomial order r ≥ 1.
    pub r: usize,
    /// 1D lattice size m·r + 1.
    pub n1: usize,
    /// 1D element mass matrix ((r+1)²) for h = 1/m.
    pub mass_e: Vec<f64>,
    /// 1D element stiffness matrix for h = 1/m.
    pub stiff_e: Vec<f64>,
}

impl TensorSpace {
    /// Build the space (uniform grid, h = 1/m per axis).
    ///
    /// # Panics
    /// If `m == 0` or `r == 0`.
    #[must_use]
    pub fn new(m: usize, r: usize) -> TensorSpace {
        assert!(m >= 1 && r >= 1, "TensorSpace needs m >= 1, r >= 1");
        let h = 1.0 / m as f64;
        let (mass_e, stiff_e) = element_matrices(r, h);
        TensorSpace {
            m,
            r,
            n1: m * r + 1,
            mass_e,
            stiff_e,
        }
    }

    /// Total dof count (m·r + 1)³.
    #[must_use]
    pub fn ndof(&self) -> usize {
        self.n1 * self.n1 * self.n1
    }

    /// Global 1D lattice index of local dof `l` (Lobatto order:
    /// 0 = left vertex, 1 = right vertex, k ≥ 2 = bubbles) in cell c.
    #[must_use]
    pub fn lat1(&self, c: usize, l: usize) -> usize {
        match l {
            0 => c * self.r,
            1 => (c + 1) * self.r,
            k => c * self.r + k - 1,
        }
    }

    /// Global dof id of lattice point (i, j, k).
    #[must_use]
    pub fn gid(&self, i: usize, j: usize, k: usize) -> usize {
        (i * self.n1 + j) * self.n1 + k
    }

    /// True when lattice index `i` carries nonzero trace on a domain
    /// face along its axis (only the two endpoint vertex-chain dofs).
    #[must_use]
    pub fn on_axis_boundary(&self, i: usize) -> bool {
        i == 0 || i == self.m * self.r
    }

    /// SUM-FACTORIZED Poisson apply y = K·u (full space, no BC).
    /// Per element: y_loc = (K⊗M⊗M + M⊗K⊗M + M⊗M⊗K)·u_loc via three
    /// axis contractions per term.
    #[must_use]
    pub fn apply_stiffness(&self, u: &[f64]) -> Vec<f64> {
        use super::fma::apply_mono_dispatch;
        assert_eq!(u.len(), self.ndof(), "apply_stiffness length mismatch");
        let mut y = vec![0.0f64; self.ndof()];
        // Dispatch ONCE per apply to a const-P element loop for the
        // practical degrees: compile-time trip counts let the whole
        // gather → 3-term contraction → scatter pipeline inline,
        // unroll, and NEON/SSE-vectorize. Every monomorphized kernel
        // preserves the generic path's per-output operation ORDER
        // exactly (ascending-l accumulation, same `mul_add` operands;
        // the dropped zero-skip is semantics-free — see the kernel
        // comment), so the output is BIT-IDENTICAL and the sf-kron
        // golden does not move.
        // Routed through the fma capsule (bead a55x): on x86 with
        // avx2+fma the SAME body compiles with native fused ops
        // (baseline x86 lowers mul_add to a per-element libm CALL —
        // measured 0.026 attainment); elsewhere identical codegen to
        // calling the body directly. Bit-identical either way.
        match self.r + 1 {
            2 => apply_mono_dispatch::<2>(self, u, &mut y),
            3 => apply_mono_dispatch::<3>(self, u, &mut y),
            4 => apply_mono_dispatch::<4>(self, u, &mut y),
            5 => apply_mono_dispatch::<5>(self, u, &mut y),
            6 => apply_mono_dispatch::<6>(self, u, &mut y),
            7 => apply_mono_dispatch::<7>(self, u, &mut y),
            8 => apply_mono_dispatch::<8>(self, u, &mut y),
            _ => self.apply_gen(u, &mut y),
        }
        y
    }

    /// The const-P element loop (P = r+1 ≤ 8): scratch on the stack,
    /// contraction kernels inlined, gather/scatter as hoisted address
    /// arithmetic. Same element visitation and accumulation order as
    /// [`Self::apply_gen`] — bit-identical output.
    #[allow(clippy::inline_always)] // must inline INTO the fma capsule's target_feature fn
    #[inline(always)]
    pub(crate) fn apply_mono_body<const P: usize>(&self, u: &[f64], y: &mut [f64]) {
        const { assert!(P <= 8) }
        let n1 = self.n1;
        let (mut local, mut t1, mut t2, mut acc) =
            ([0.0f64; 512], [0.0f64; 512], [0.0f64; 512], [0.0f64; 512]);
        let (local, t1, t2) = (
            &mut local[..P * P * P],
            &mut t1[..P * P * P],
            &mut t2[..P * P * P],
        );
        let acc = &mut acc[..P * P * P];
        let (mut gx, mut gy, mut gz) = ([0usize; 8], [0usize; 8], [0usize; 8]);
        for cx in 0..self.m {
            for (l, g) in gx[..P].iter_mut().enumerate() {
                *g = self.lat1(cx, l);
            }
            for cy in 0..self.m {
                for (l, g) in gy[..P].iter_mut().enumerate() {
                    *g = self.lat1(cy, l);
                }
                for cz in 0..self.m {
                    for (l, g) in gz[..P].iter_mut().enumerate() {
                        *g = self.lat1(cz, l);
                    }
                    // Gather.
                    for lx in 0..P {
                        for ly in 0..P {
                            let base = (gx[lx] * n1 + gy[ly]) * n1;
                            let row = &mut local[(lx * P + ly) * P..(lx * P + ly + 1) * P];
                            for (lz, v) in row.iter_mut().enumerate() {
                                *v = u[base + gz[lz]];
                            }
                        }
                    }
                    // Three Kronecker terms; contract_*_p applies a
                    // dense (P×P) 1D matrix along one axis.
                    acc.fill(0.0);
                    for term in 0..3 {
                        let (ax, ay, az) = match term {
                            0 => (&self.stiff_e, &self.mass_e, &self.mass_e),
                            1 => (&self.mass_e, &self.stiff_e, &self.mass_e),
                            _ => (&self.mass_e, &self.mass_e, &self.stiff_e),
                        };
                        contract_x_p::<P>(ax, local, t1);
                        contract_y_p::<P>(ay, t1, t2);
                        contract_z_p::<P>(az, t2, t1);
                        for (a, t) in acc.iter_mut().zip(&*t1) {
                            *a += *t;
                        }
                    }
                    // Scatter-add (fixed element order).
                    for lx in 0..P {
                        for ly in 0..P {
                            let base = (gx[lx] * n1 + gy[ly]) * n1;
                            let row = &acc[(lx * P + ly) * P..(lx * P + ly + 1) * P];
                            for (lz, &v) in row.iter().enumerate() {
                                y[base + gz[lz]] += v;
                            }
                        }
                    }
                }
            }
        }
    }

    /// The runtime-p fallback (r ≥ 8) — the original loop structure.
    fn apply_gen(&self, u: &[f64], y: &mut [f64]) {
        let p = self.r + 1;
        let mut local = vec![0.0f64; p * p * p];
        let mut t1 = vec![0.0f64; p * p * p];
        let mut t2 = vec![0.0f64; p * p * p];
        let mut acc = vec![0.0f64; p * p * p];
        for cx in 0..self.m {
            for cy in 0..self.m {
                for cz in 0..self.m {
                    // Gather.
                    for lx in 0..p {
                        let gi = self.lat1(cx, lx);
                        for ly in 0..p {
                            let gj = self.lat1(cy, ly);
                            for lz in 0..p {
                                let gk = self.lat1(cz, lz);
                                local[(lx * p + ly) * p + lz] = u[self.gid(gi, gj, gk)];
                            }
                        }
                    }
                    acc.fill(0.0);
                    for term in 0..3 {
                        let (ax, ay, az) = match term {
                            0 => (&self.stiff_e, &self.mass_e, &self.mass_e),
                            1 => (&self.mass_e, &self.stiff_e, &self.mass_e),
                            _ => (&self.mass_e, &self.mass_e, &self.stiff_e),
                        };
                        contract_x_gen(ax, &local, &mut t1, p);
                        contract_y_gen(ay, &t1, &mut t2, p);
                        contract_z_gen(az, &t2, &mut t1, p);
                        for (a, t) in acc.iter_mut().zip(&t1) {
                            *a += t;
                        }
                    }
                    // Scatter-add (fixed element order).
                    for lx in 0..p {
                        let gi = self.lat1(cx, lx);
                        for ly in 0..p {
                            let gj = self.lat1(cy, ly);
                            for lz in 0..p {
                                let gk = self.lat1(cz, lz);
                                y[self.gid(gi, gj, gk)] += acc[(lx * p + ly) * p + lz];
                            }
                        }
                    }
                }
            }
        }
    }

    /// Assembled 1D global matrices (mass, stiffness) — dense n₁×n₁,
    /// the reference the Kronecker comparison and the Jacobi diagonal
    /// are built from.
    #[must_use]
    pub fn assembled_1d(&self) -> (Vec<f64>, Vec<f64>) {
        let p = self.r + 1;
        let n1 = self.n1;
        let mut mass = vec![0.0f64; n1 * n1];
        let mut stiff = vec![0.0f64; n1 * n1];
        for c in 0..self.m {
            for li in 0..p {
                let gi = self.lat1(c, li);
                for lj in 0..p {
                    let gj = self.lat1(c, lj);
                    mass[gi * n1 + gj] += self.mass_e[li * p + lj];
                    stiff[gi * n1 + gj] += self.stiff_e[li * p + lj];
                }
            }
        }
        (mass, stiff)
    }

    /// Exact operator diagonal from the Kronecker structure:
    /// diag[(i,j,k)] = Kd·Md·Md + Md·Kd·Md + Md·Md·Kd (the matrix-free
    /// Jacobi preconditioner — P6: never assemble what we can apply).
    #[must_use]
    pub fn stiffness_diagonal(&self) -> Vec<f64> {
        let (mass, stiff) = self.assembled_1d();
        let n1 = self.n1;
        let md: Vec<f64> = (0..n1).map(|i| mass[i * n1 + i]).collect();
        let kd: Vec<f64> = (0..n1).map(|i| stiff[i * n1 + i]).collect();
        let mut diag = vec![0.0f64; self.ndof()];
        for i in 0..n1 {
            for j in 0..n1 {
                for k in 0..n1 {
                    diag[self.gid(i, j, k)] = kd[i].mul_add(
                        md[j] * md[k],
                        md[i].mul_add(kd[j] * md[k], md[i] * md[j] * kd[k]),
                    );
                }
            }
        }
        diag
    }

    /// Load vector b_i = ∫ f·φ_i by per-element tensor quadrature
    /// (r+4 points per axis: quadrature error far below discretization
    /// error for smooth f).
    #[must_use]
    pub fn load<F: Fn([f64; 3]) -> f64>(&self, f: &F) -> Vec<f64> {
        let p = self.r + 1;
        let h = 1.0 / self.m as f64;
        let nq = self.r + 4;
        let (qx, qw) = gauss_legendre(nq);
        // Basis values at quadrature points (shared across elements).
        let shapes: Vec<Vec<f64>> = qx.iter().map(|&x| lobatto_shapes(self.r, x).0).collect();
        let mut b = vec![0.0f64; self.ndof()];
        let jac = (h / 2.0) * (h / 2.0) * (h / 2.0);
        for cx in 0..self.m {
            for cy in 0..self.m {
                for cz in 0..self.m {
                    for (qi, &xq) in qx.iter().enumerate() {
                        let px = (cx as f64).mul_add(h, (xq + 1.0) * h / 2.0);
                        for (qj, &yq) in qx.iter().enumerate() {
                            let py = (cy as f64).mul_add(h, (yq + 1.0) * h / 2.0);
                            for (qk, &zq) in qx.iter().enumerate() {
                                let pz = (cz as f64).mul_add(h, (zq + 1.0) * h / 2.0);
                                let w = qw[qi] * qw[qj] * qw[qk] * jac * f([px, py, pz]);
                                for lx in 0..p {
                                    let gi = self.lat1(cx, lx);
                                    let sx = shapes[qi][lx];
                                    for ly in 0..p {
                                        let gj = self.lat1(cy, ly);
                                        let sxy = sx * shapes[qj][ly];
                                        for (lz, sz) in shapes[qk].iter().enumerate() {
                                            let gk = self.lat1(cz, lz);
                                            b[self.gid(gi, gj, gk)] += w * sxy * sz;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        b
    }

    /// L2 error ‖u_h − u‖ by per-element tensor quadrature.
    #[must_use]
    pub fn l2_error<F: Fn([f64; 3]) -> f64>(&self, u_dofs: &[f64], u_exact: &F) -> f64 {
        let p = self.r + 1;
        let h = 1.0 / self.m as f64;
        let nq = self.r + 4;
        let (qx, qw) = gauss_legendre(nq);
        let shapes: Vec<Vec<f64>> = qx.iter().map(|&x| lobatto_shapes(self.r, x).0).collect();
        let jac = (h / 2.0) * (h / 2.0) * (h / 2.0);
        let mut total = 0.0f64;
        for cx in 0..self.m {
            for cy in 0..self.m {
                for cz in 0..self.m {
                    for (qi, &xq) in qx.iter().enumerate() {
                        let px = (cx as f64).mul_add(h, (xq + 1.0) * h / 2.0);
                        for (qj, &yq) in qx.iter().enumerate() {
                            let py = (cy as f64).mul_add(h, (yq + 1.0) * h / 2.0);
                            for (qk, &zq) in qx.iter().enumerate() {
                                let pz = (cz as f64).mul_add(h, (zq + 1.0) * h / 2.0);
                                let mut uh = 0.0f64;
                                for lx in 0..p {
                                    let gi = self.lat1(cx, lx);
                                    for ly in 0..p {
                                        let gj = self.lat1(cy, ly);
                                        for lz in 0..p {
                                            let gk = self.lat1(cz, lz);
                                            uh += u_dofs[self.gid(gi, gj, gk)]
                                                * shapes[qi][lx]
                                                * shapes[qj][ly]
                                                * shapes[qk][lz];
                                        }
                                    }
                                }
                                let e = uh - u_exact([px, py, pz]);
                                total += qw[qi] * qw[qj] * qw[qk] * jac * e * e;
                            }
                        }
                    }
                }
            }
        }
        fs_math::det::sqrt(total)
    }

    /// Interior-dof mask (homogeneous Dirichlet on all six faces): a
    /// 3D basis function has nonzero boundary trace iff any axis
    /// lattice index is an endpoint vertex dof.
    #[must_use]
    pub fn interior_mask(&self) -> Vec<bool> {
        let n1 = self.n1;
        let mut mask = vec![false; self.ndof()];
        for i in 0..n1 {
            for j in 0..n1 {
                for k in 0..n1 {
                    mask[self.gid(i, j, k)] = !(self.on_axis_boundary(i)
                        || self.on_axis_boundary(j)
                        || self.on_axis_boundary(k));
                }
            }
        }
        mask
    }
}

// CONST-P monomorphized contraction kernels for the practical degrees
// (p = r+1 in 2..=8): compile-time trip counts let LLVM fully unroll
// and NEON/SSE-vectorize the stride-1 inner loops, which the runtime-p
// `_gen` fallbacks never achieve. Every kernel preserves the
// fallback's per-output accumulation ORDER exactly: ascending l,
// row-major positions, same `mul_add` operands. The fallbacks' ail/ajl
// zero-skip is DROPPED here (branches block straight-line register
// allocation of the unrolled GEMMs); that is still bit-identical:
// executing the skipped op computes fma(±0·s, acc) = acc, because acc
// is never -0.0 (it starts +0.0, and round-to-nearest addition onto
// ±0/nonzero never yields -0.0 from a +0.0 start).
//
// `inline(always)`: these are the innermost hot kernels of the whole
// apply and MUST fuse into the monomorphized element loop — measured
// on the perf lane, not assumed.
// REGISTER-ARRAY ACCUMULATORS (bead cwjn, the a55x handoff's lead): each
// output row accumulates in a local `[f64; P]` and stores ONCE, instead
// of fma-ing through `&mut dst` memory. Bit-identical per element — the
// chain is still "+0.0, then ascending-l mul_add" — but the value-array
// form removes the load/store round-trip and every aliasing question
// between passes, which is exactly what x86 SLP needs to pack the
// unrolled const-P loops into ymm ops (the a55x objdump found all 69
// vfmadds SCALAR on x86 in the through-memory form, while aarch64
// autovectorized it). The `super::fma` target_feature capsule one level
// up keeps mul_add an inline vfmadd on x86.
//
// `inline(always)`: these are the innermost hot kernels of the whole
// apply and MUST fuse into the monomorphized element loop — measured
// on the perf lane, not assumed.
#[allow(clippy::inline_always)]
#[inline(always)]
fn contract_x_p<const P: usize>(a: &[f64], src: &[f64], dst: &mut [f64]) {
    for i in 0..P {
        let a_row = &a[i * P..(i + 1) * P];
        for j in 0..P {
            let mut acc = [0.0f64; P];
            for l in 0..P {
                let ail = a_row[l];
                let srow = &src[(l * P + j) * P..(l * P + j + 1) * P];
                for k in 0..P {
                    acc[k] = ail.mul_add(srow[k], acc[k]);
                }
            }
            dst[(i * P + j) * P..(i * P + j + 1) * P].copy_from_slice(&acc);
        }
    }
}

#[allow(clippy::inline_always)]
#[inline(always)]
fn contract_y_p<const P: usize>(a: &[f64], src: &[f64], dst: &mut [f64]) {
    for i in 0..P {
        let sblock = &src[i * P * P..(i + 1) * P * P];
        for j in 0..P {
            let a_row = &a[j * P..(j + 1) * P];
            let mut acc = [0.0f64; P];
            for l in 0..P {
                let ajl = a_row[l];
                let srow = &sblock[l * P..(l + 1) * P];
                for k in 0..P {
                    acc[k] = ajl.mul_add(srow[k], acc[k]);
                }
            }
            dst[(i * P + j) * P..(i * P + j + 1) * P].copy_from_slice(&acc);
        }
    }
}

#[allow(clippy::inline_always)]
#[inline(always)]
fn contract_z_p<const P: usize>(a: &[f64], src: &[f64], dst: &mut [f64]) {
    // Transpose the 1D matrix (values unchanged) so the k loop is
    // stride-1 in both the matrix and the accumulator; each output
    // still starts at 0.0 and accumulates over ascending l — the same
    // op sequence per element as the k-inner fallback, bit-identical.
    let mut at = [0.0f64; 64]; // P <= 8, so P*P <= 64
    for k in 0..P {
        for l in 0..P {
            at[l * P + k] = a[k * P + l];
        }
    }
    for (drow, srow) in dst
        .as_chunks_mut::<P>()
        .0
        .iter_mut()
        .zip(src.as_chunks::<P>().0)
    {
        let mut acc = [0.0f64; P];
        for l in 0..P {
            let sl = srow[l];
            let arow = &at[l * P..(l + 1) * P];
            for k in 0..P {
                acc[k] = arow[k].mul_add(sl, acc[k]);
            }
        }
        *drow = acc;
    }
}

fn contract_x_gen(a: &[f64], src: &[f64], dst: &mut [f64], p: usize) {
    dst.fill(0.0);
    for i in 0..p {
        for l in 0..p {
            let ail = a[i * p + l];
            if ail != 0.0 {
                for j in 0..p {
                    for k in 0..p {
                        dst[(i * p + j) * p + k] =
                            ail.mul_add(src[(l * p + j) * p + k], dst[(i * p + j) * p + k]);
                    }
                }
            }
        }
    }
}

fn contract_y_gen(a: &[f64], src: &[f64], dst: &mut [f64], p: usize) {
    dst.fill(0.0);
    for i in 0..p {
        for j in 0..p {
            for l in 0..p {
                let ajl = a[j * p + l];
                if ajl != 0.0 {
                    for k in 0..p {
                        dst[(i * p + j) * p + k] =
                            ajl.mul_add(src[(i * p + l) * p + k], dst[(i * p + j) * p + k]);
                    }
                }
            }
        }
    }
}

fn contract_z_gen(a: &[f64], src: &[f64], dst: &mut [f64], p: usize) {
    dst.fill(0.0);
    for i in 0..p {
        for j in 0..p {
            for k in 0..p {
                let mut acc = 0.0f64;
                for l in 0..p {
                    acc = a[k * p + l].mul_add(src[(i * p + j) * p + l], acc);
                }
                dst[(i * p + j) * p + k] = acc;
            }
        }
    }
}

/// Matrix-free Jacobi-PCG on the interior dofs of a homogeneous
/// Dirichlet problem: solve K u = b with `apply` the FULL-space
/// operator (boundary dofs zeroed each application). Returns
/// (iterations, converged).
pub fn pcg_matfree<A: Fn(&[f64]) -> Vec<f64>>(
    apply: &A,
    b: &[f64],
    x: &mut [f64],
    mask: &[bool],
    diag: &[f64],
    tol: f64,
    max_iters: usize,
) -> (usize, bool) {
    let n = b.len();
    let project = |v: &mut [f64]| {
        for (vi, &m) in v.iter_mut().zip(mask) {
            if !m {
                *vi = 0.0;
            }
        }
    };
    let bnorm = fs_math::det::sqrt(
        b.iter()
            .zip(mask)
            .filter(|(_, m)| **m)
            .map(|(v, _)| v * v)
            .sum::<f64>(),
    )
    .max(f64::MIN_POSITIVE);
    let mut r: Vec<f64> = {
        let mut ax = apply(x);
        project(&mut ax);
        b.iter().zip(&ax).map(|(bi, ai)| bi - ai).collect()
    };
    project(&mut r);
    let mut z: Vec<f64> = r.iter().zip(diag).map(|(ri, di)| ri / di).collect();
    project(&mut z);
    let mut p_dir = z.clone();
    let mut rz: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
    for it in 0..max_iters {
        let rnorm = fs_math::det::sqrt(r.iter().map(|v| v * v).sum::<f64>());
        if rnorm / bnorm < tol {
            return (it, true);
        }
        let mut ap = apply(&p_dir);
        project(&mut ap);
        let pap: f64 = p_dir.iter().zip(&ap).map(|(a, b)| a * b).sum();
        let alpha = rz / pap;
        for i in 0..n {
            x[i] = alpha.mul_add(p_dir[i], x[i]);
            r[i] = alpha.mul_add(-ap[i], r[i]);
        }
        for i in 0..n {
            z[i] = if mask[i] { r[i] / diag[i] } else { 0.0 };
        }
        let rz_new: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
        let beta = rz_new / rz;
        rz = rz_new;
        for i in 0..n {
            p_dir[i] = beta.mul_add(p_dir[i], z[i]);
        }
    }
    (max_iters, false)
}

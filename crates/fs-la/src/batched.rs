//! Batched small dense LA (plan §6.1) — the ACTUAL hot loop of
//! FEM-adjacent work: element matrices are 6×6…48×48 and dominate
//! assembly-adjacent compute.
//!
//! Layout doctrine: interleave ACROSS the batch. [`BatchMat`] stores
//! one 128-byte-aligned SoA "plane" per matrix entry (i, j), each
//! plane holding that entry for every matrix in the batch — SIMD lanes
//! run across ELEMENTS, never within one tiny matrix. Because each
//! matrix's arithmetic is independent and every loop is fixed-order
//! `mul_add`, results are bit-deterministic by construction AND
//! independent of batch size and position (tested: a matrix computed
//! in a batch of N is bitwise-equal to the same matrix in a batch of
//! 1). Storage is `fs_soa::FieldBuf` — the wf9.5 substrate.
//!
//! Per-element failures (non-SPD pivot, exactly singular column) are
//! FLAGGED and the batch CONTINUES: flagged matrices' outputs are
//! unspecified-but-finite, the flag list is authoritative. Error
//! taxonomy mirrors `factor::FactorError`.

use crate::factor::FactorError;
use fs_math::det;
use fs_soa::FieldBuf;

/// Plane stride granularity: 16 f64 = 128 bytes, so every plane start
/// stays 128-byte aligned given the buffer base is.
const STRIDE_QUANTUM: usize = 16;

const fn stride_for(n: usize) -> usize {
    let s = n.div_ceil(STRIDE_QUANTUM) * STRIDE_QUANTUM;
    // Break power-of-two plane spacing (bead 9ekv): 2ᵖ strides put
    // every plane in the same L1/L2 cache sets and the multi-stream
    // tile kernels thrash — measured 2× slower at k ∈ {8, 16, 32}
    // batch sizes than at non-power-of-two neighbours. One extra
    // quantum of padding is pure LAYOUT: entry values and operation
    // order are untouched (bit-neutral, covered by the membership-
    // invariance and golden tests).
    if s >= 4096 && s.is_power_of_two() {
        s + STRIDE_QUANTUM
    } else {
        s
    }
}

fn make_buf(planes: usize, stride: usize) -> FieldBuf<f64> {
    let total = planes * stride;
    let mut data = FieldBuf::with_capacity(total);
    for _ in 0..total {
        data.push(0.0);
    }
    data
}

/// A batch of `n` dense k×k matrices in entry-plane SoA layout: plane
/// p = i·k + j holds entry (i, j) of every matrix, 128-byte aligned.
#[derive(Debug, Clone)]
pub struct BatchMat {
    data: FieldBuf<f64>,
    k: usize,
    n: usize,
    stride: usize,
}

impl BatchMat {
    /// Batch of `n` zero k×k matrices.
    ///
    /// # Panics
    /// If `k == 0` or `n == 0`.
    #[must_use]
    pub fn zeros(k: usize, n: usize) -> BatchMat {
        assert!(k > 0 && n > 0, "empty batch shapes are programmer errors");
        let stride = stride_for(n);
        BatchMat {
            data: make_buf(k * k, stride),
            k,
            n,
            stride,
        }
    }

    /// Batch filled from `f(m, i, j)` (matrix index, row, column).
    #[must_use]
    pub fn from_fn<F: FnMut(usize, usize, usize) -> f64>(k: usize, n: usize, mut f: F) -> BatchMat {
        let mut out = BatchMat::zeros(k, n);
        for i in 0..k {
            for j in 0..k {
                let plane = out.plane_mut(i, j);
                for (m, slot) in plane.iter_mut().enumerate() {
                    *slot = f(m, i, j);
                }
            }
        }
        out
    }

    /// Matrix dimension k.
    #[must_use]
    pub const fn k(&self) -> usize {
        self.k
    }

    /// Batch size n.
    #[must_use]
    pub const fn batch_len(&self) -> usize {
        self.n
    }

    /// Entry (i, j) of matrix `m`.
    #[must_use]
    pub fn get(&self, m: usize, i: usize, j: usize) -> f64 {
        self.plane(i, j)[m]
    }

    /// Set entry (i, j) of matrix `m`.
    pub fn set(&mut self, m: usize, i: usize, j: usize, v: f64) {
        self.plane_mut(i, j)[m] = v;
    }

    /// The (i, j) entry plane across the whole batch (len = batch).
    #[must_use]
    pub fn plane(&self, i: usize, j: usize) -> &[f64] {
        assert!(i < self.k && j < self.k, "plane index out of bounds");
        let p = i * self.k + j;
        &self.data.as_slice()[p * self.stride..p * self.stride + self.n]
    }

    /// Mutable (i, j) entry plane across the whole batch.
    pub fn plane_mut(&mut self, i: usize, j: usize) -> &mut [f64] {
        assert!(i < self.k && j < self.k, "plane index out of bounds");
        let p = i * self.k + j;
        let (stride, n) = (self.stride, self.n);
        &mut self.data.as_mut_slice()[p * stride..p * stride + n]
    }

    /// Two distinct mutable planes at once (target, source pattern).
    ///
    /// # Panics
    /// If the plane indices coincide.
    pub fn planes_mut2(
        &mut self,
        a: (usize, usize),
        b: (usize, usize),
    ) -> (&mut [f64], &mut [f64]) {
        let pa = a.0 * self.k + a.1;
        let pb = b.0 * self.k + b.1;
        assert!(pa != pb, "planes_mut2 requires distinct planes");
        let (stride, n) = (self.stride, self.n);
        let s = self.data.as_mut_slice();
        if pa < pb {
            let (lo, hi) = s.split_at_mut(pb * stride);
            (&mut lo[pa * stride..pa * stride + n], &mut hi[..n])
        } else {
            let (lo, hi) = s.split_at_mut(pa * stride);
            (&mut hi[..n], &mut lo[pb * stride..pb * stride + n])
        }
    }

    /// Gather matrix `m` into row-major AoS form (k×k).
    #[must_use]
    pub fn gather(&self, m: usize) -> Vec<f64> {
        let k = self.k;
        let mut out = vec![0.0f64; k * k];
        for i in 0..k {
            for j in 0..k {
                out[i * k + j] = self.get(m, i, j);
            }
        }
        out
    }

    /// Scatter a row-major AoS matrix into slot `m`.
    pub fn scatter(&mut self, m: usize, a: &[f64]) {
        let k = self.k;
        assert_eq!(a.len(), k * k, "scatter shape mismatch");
        for i in 0..k {
            for j in 0..k {
                self.set(m, i, j, a[i * k + j]);
            }
        }
    }
}

/// A batch of `n` k-vectors in component-plane SoA layout.
#[derive(Debug, Clone)]
pub struct BatchVec {
    data: FieldBuf<f64>,
    k: usize,
    n: usize,
    stride: usize,
}

impl BatchVec {
    /// Batch of `n` zero k-vectors.
    ///
    /// # Panics
    /// If `k == 0` or `n == 0`.
    #[must_use]
    pub fn zeros(k: usize, n: usize) -> BatchVec {
        assert!(k > 0 && n > 0, "empty batch shapes are programmer errors");
        let stride = stride_for(n);
        BatchVec {
            data: make_buf(k, stride),
            k,
            n,
            stride,
        }
    }

    /// Batch filled from `f(m, i)`.
    #[must_use]
    pub fn from_fn<F: FnMut(usize, usize) -> f64>(k: usize, n: usize, mut f: F) -> BatchVec {
        let mut out = BatchVec::zeros(k, n);
        for i in 0..k {
            let plane = out.plane_mut(i);
            for (m, slot) in plane.iter_mut().enumerate() {
                *slot = f(m, i);
            }
        }
        out
    }

    /// Vector dimension k.
    #[must_use]
    pub const fn k(&self) -> usize {
        self.k
    }

    /// Batch size n.
    #[must_use]
    pub const fn batch_len(&self) -> usize {
        self.n
    }

    /// Component i of vector `m`.
    #[must_use]
    pub fn get(&self, m: usize, i: usize) -> f64 {
        self.plane(i)[m]
    }

    /// Set component i of vector `m`.
    pub fn set(&mut self, m: usize, i: usize, v: f64) {
        self.plane_mut(i)[m] = v;
    }

    /// Component-i plane across the batch.
    #[must_use]
    pub fn plane(&self, i: usize) -> &[f64] {
        assert!(i < self.k, "plane index out of bounds");
        &self.data.as_slice()[i * self.stride..i * self.stride + self.n]
    }

    /// Mutable component-i plane across the batch.
    pub fn plane_mut(&mut self, i: usize) -> &mut [f64] {
        assert!(i < self.k, "plane index out of bounds");
        let (stride, n) = (self.stride, self.n);
        &mut self.data.as_mut_slice()[i * stride..i * stride + n]
    }

    /// Gather vector `m` into AoS form.
    #[must_use]
    pub fn gather(&self, m: usize) -> Vec<f64> {
        (0..self.k).map(|i| self.get(m, i)).collect()
    }
}

// -------------------------------------------------------------------- GEMM

/// Batched GEMM: per matrix, C ← α·A·B + β·C (β = 0 overwrites C —
/// the BLAS uninitialized-output convention, matching `gemm::gemm_f64`).
/// Fixed l-ascending accumulation with `mul_add`: bit-deterministic
/// and batch-size independent. Size classes {4, 6, 8, 12, 16, 24, 32,
/// 48} dispatch to monomorphized kernels; other k take the same code
/// generically.
///
/// # Panics
/// On shape or batch-length mismatch.
pub fn batch_gemm(alpha: f64, a: &BatchMat, b: &BatchMat, beta: f64, c: &mut BatchMat) {
    let k = a.k;
    assert!(b.k == k && c.k == k, "batch_gemm dimension mismatch");
    assert!(a.n == b.n && a.n == c.n, "batch_gemm batch-length mismatch");
    match k {
        4 => gemm_sized::<4>(alpha, a, b, beta, c),
        6 => gemm_sized::<6>(alpha, a, b, beta, c),
        8 => gemm_sized::<8>(alpha, a, b, beta, c),
        12 => gemm_sized::<12>(alpha, a, b, beta, c),
        16 => gemm_sized::<16>(alpha, a, b, beta, c),
        24 => gemm_sized::<24>(alpha, a, b, beta, c),
        32 => gemm_sized::<32>(alpha, a, b, beta, c),
        48 => gemm_sized::<48>(alpha, a, b, beta, c),
        _ => gemm_generic(alpha, a, b, beta, c, k),
    }
}

fn gemm_sized<const K: usize>(alpha: f64, a: &BatchMat, b: &BatchMat, beta: f64, c: &mut BatchMat) {
    gemm_generic(alpha, a, b, beta, c, K);
}

/// Lanes per batch chunk (bead 9ekv): the accumulator chunk stays
/// L1-resident and the active plane slices L2-resident instead of
/// streaming the whole batch from memory k² times. PURE lane
/// partitioning — per-element operation order is untouched, so this is
/// bit-neutral (the membership-invariance contract, tested).
///
/// TUNED EMPIRICALLY on M4 Pro (release, median-of-9, batched_perf):
/// 256 → 7–18% of peak (the 16·256·8 B = 32 KB accumulator tile plus
/// streamed planes overran L1 residency); 128 worse; 64 → 9–23%;
/// 32 → 9–26% (best all-rounder); 16 ≈ 32 but weaker at k = 24/48.
/// The remaining gap to the 60% band lives in the btile kernel's
/// inner loop (the fs-simd capsule lane), not in this partitioning.
const MBLK: usize = 32;

fn gemm_generic(alpha: f64, a: &BatchMat, b: &BatchMat, beta: f64, c: &mut BatchMat, k: usize) {
    let n = a.n;
    let btile = fs_simd::ops().btile4x4p_f64;
    let mut tile = vec![0.0f64; 16 * MBLK];
    let mut acc = [0.0f64; MBLK];
    // PACKED operands (bead 9ekv slice 2): per chunk, copy the lane
    // block of every plane into l-CONTIGUOUS scratch — A i-major
    // (a_pack[(i·k + l)·mb]), B j-major (b_pack[(j·k + l)·mb]) — so the
    // tile kernel's walks stream sequentially instead of jumping
    // k·stride pages (measured TLB/latency bound). Packing is 2k²·mb
    // copies against 2k³·mb flops (1/k overhead) and is PURE DATA
    // MOVEMENT: per-element operation order is unchanged, so outputs
    // stay bitwise-identical (membership invariance + golden, tested).
    let mut a_pack = vec![0.0f64; k * k * MBLK];
    let mut b_pack = vec![0.0f64; k * k * MBLK];
    let kt = k - k % 4; // tile-covered leading square
    let mut m0 = 0;
    while m0 < n {
        let mb = MBLK.min(n - m0);
        for i in 0..k {
            for l in 0..k {
                a_pack[(i * k + l) * mb..(i * k + l) * mb + mb]
                    .copy_from_slice(&a.plane(i, l)[m0..m0 + mb]);
            }
        }
        for j in 0..k {
            for l in 0..k {
                b_pack[(j * k + l) * mb..(j * k + l) * mb + mb]
                    .copy_from_slice(&b.plane(l, j)[m0..m0 + mb]);
            }
        }
        for i0 in (0..kt).step_by(4) {
            for j0 in (0..kt).step_by(4) {
                let td = &mut tile[..16 * mb];
                btile(&a_pack, &b_pack, i0, j0, k, mb, td);
                for ti in 0..4 {
                    for tj in 0..4 {
                        let trow = &td[(ti * 4 + tj) * mb..(ti * 4 + tj + 1) * mb];
                        let cp = &mut c.plane_mut(i0 + ti, j0 + tj)[m0..m0 + mb];
                        write_back(alpha, beta, trow, cp);
                    }
                }
            }
        }
        // Tails (k % 4 rows and columns): the plane-at-a-time loop on
        // the ORIGINAL planes (identical ops — bitwise-consistent).
        for i in 0..k {
            let j_start = if i < kt { kt } else { 0 };
            for j in j_start..k {
                let acc = &mut acc[..mb];
                acc.fill(0.0);
                for l in 0..k {
                    let ap = &a.plane(i, l)[m0..m0 + mb];
                    let bp = &b.plane(l, j)[m0..m0 + mb];
                    for ((s, &am), &bm) in acc.iter_mut().zip(ap).zip(bp) {
                        *s = am.mul_add(bm, *s);
                    }
                }
                let cp = &mut c.plane_mut(i, j)[m0..m0 + mb];
                write_back(alpha, beta, acc, cp);
            }
        }
        m0 += MBLK;
    }
}

/// α/β application per output plane chunk (β = 0 overwrites — the BLAS
/// convention; identical op per element on both the tile and tail
/// paths, so the two paths are bitwise-consistent).
fn write_back(alpha: f64, beta: f64, acc: &[f64], cp: &mut [f64]) {
    if beta == 0.0 {
        for (cm, &s) in cp.iter_mut().zip(acc) {
            *cm = alpha * s;
        }
    } else {
        for (cm, &s) in cp.iter_mut().zip(acc) {
            *cm = alpha.mul_add(s, beta * *cm);
        }
    }
}

// ------------------------------------------------------- det / inverse ≤ 4

/// Batched determinants, closed form for k ∈ {1, 2, 3, 4} (Jacobian
/// hot path). Larger k: factor with [`batch_lu`] instead.
///
/// # Panics
/// If `k > 4`.
#[must_use]
pub fn batch_det(a: &BatchMat) -> Vec<f64> {
    let n = a.n;
    let mut out = vec![0.0f64; n];
    match a.k {
        1 => out.copy_from_slice(a.plane(0, 0)),
        2 => {
            let (a00, a01, a10, a11) = (a.plane(0, 0), a.plane(0, 1), a.plane(1, 0), a.plane(1, 1));
            for (m, o) in out.iter_mut().enumerate() {
                *o = a00[m].mul_add(a11[m], -(a01[m] * a10[m]));
            }
        }
        3 => {
            for (m, o) in out.iter_mut().enumerate() {
                *o = det3(&gather3(a, m));
            }
        }
        4 => {
            for (m, o) in out.iter_mut().enumerate() {
                *o = det4(&a.gather(m));
            }
        }
        _ => panic!("batch_det closed form covers k <= 4; use batch_lu for larger"),
    }
    out
}

fn gather3(a: &BatchMat, m: usize) -> [f64; 9] {
    let mut g = [0.0f64; 9];
    for i in 0..3 {
        for j in 0..3 {
            g[i * 3 + j] = a.get(m, i, j);
        }
    }
    g
}

fn det3(g: &[f64; 9]) -> f64 {
    let c0 = g[4].mul_add(g[8], -(g[5] * g[7]));
    let c1 = g[3].mul_add(g[8], -(g[5] * g[6]));
    let c2 = g[3].mul_add(g[7], -(g[4] * g[6]));
    g[0].mul_add(c0, g[2].mul_add(c2, -(g[1] * c1)))
}

fn det4(g: &[f64]) -> f64 {
    // Cofactor expansion along row 0 with 3×3 minors, fixed order.
    let mut total = 0.0f64;
    for j in 0..4 {
        let mut minor = [0.0f64; 9];
        let mut idx = 0;
        for r in 1..4 {
            for c in 0..4 {
                if c != j {
                    minor[idx] = g[r * 4 + c];
                    idx += 1;
                }
            }
        }
        let sign = if j % 2 == 0 { 1.0 } else { -1.0 };
        total = (sign * g[j]).mul_add(det3(&minor), total);
    }
    total
}

/// Batched inverses, closed form (adjugate/determinant) for
/// k ∈ {1, 2, 3, 4}. Exactly-zero determinants are flagged
/// `Singular` (that matrix's output is unspecified); near-singular
/// inputs honestly produce large entries — certified bounds are
/// fs-ivl territory.
///
/// # Panics
/// If `k > 4` or shapes mismatch.
#[must_use = "the flag list is the only record of singular members"]
pub fn batch_inv(a: &BatchMat, out: &mut BatchMat) -> Vec<(usize, FactorError)> {
    assert!(a.k == out.k && a.n == out.n, "batch_inv shape mismatch");
    assert!(
        a.k <= 4,
        "batch_inv closed form covers k <= 4; use batch_lu for larger"
    );
    let n = a.n;
    let mut flags = Vec::new();
    let dets = batch_det(a);
    for (m, &d) in dets.iter().enumerate().take(n) {
        if d == 0.0 {
            flags.push((m, FactorError::Singular { index: 0 }));
            continue;
        }
        let g = a.gather(m);
        let inv = inv_small(&g, a.k, d);
        out.scatter(m, &inv);
    }
    flags
}

fn inv_small(g: &[f64], k: usize, d: f64) -> Vec<f64> {
    let mut out = vec![0.0f64; k * k];
    match k {
        1 => out[0] = 1.0 / d,
        2 => {
            out[0] = g[3] / d;
            out[1] = -g[1] / d;
            out[2] = -g[2] / d;
            out[3] = g[0] / d;
        }
        _ => {
            // Adjugate via cofactors (k = 3, 4): out[j*k+i] = C_ij / d.
            for i in 0..k {
                for j in 0..k {
                    let minor = minor_of(g, k, i, j);
                    let sign = if (i + j) % 2 == 0 { 1.0 } else { -1.0 };
                    let cof = sign
                        * if k == 3 {
                            minor[0].mul_add(minor[3], -(minor[1] * minor[2]))
                        } else {
                            let mm: [f64; 9] = minor[..9].try_into().expect("3x3 minor of 4x4");
                            det3(&mm)
                        };
                    out[j * k + i] = cof / d;
                }
            }
        }
    }
    out
}

fn minor_of(g: &[f64], k: usize, i: usize, j: usize) -> Vec<f64> {
    let mut minor = Vec::with_capacity((k - 1) * (k - 1));
    for r in 0..k {
        for c in 0..k {
            if r != i && c != j {
                minor.push(g[r * k + c]);
            }
        }
    }
    minor
}

// ---------------------------------------------------------------- Cholesky

/// Batched Cholesky A = L·Lᵀ for SPD members. Non-positive pivots are
/// flagged `NotSpd { index: diagonal step }` and that matrix's factor
/// is unspecified-but-finite (the pivot is replaced by 1.0 so the
/// batch continues without NaN storms). Fixed-order `mul_add`
/// accumulation: bit-deterministic and batch-size independent.
#[must_use = "the flag list is the only record of non-SPD members"]
pub fn batch_cholesky(a: &BatchMat) -> (BatchMat, Vec<(usize, FactorError)>) {
    let (k, n) = (a.k, a.n);
    let mut l = BatchMat::zeros(k, n);
    let mut flags = Vec::new();
    let mut flagged = vec![false; n];
    for j in 0..k {
        // d[m] = A[j][j] − Σ_{p<j} L[j][p]²
        let mut d: Vec<f64> = a.plane(j, j).to_vec();
        for p in 0..j {
            let ljp = l.plane(j, p);
            for m in 0..n {
                d[m] = ljp[m].mul_add(-ljp[m], d[m]);
            }
        }
        for (m, dm) in d.iter_mut().enumerate() {
            if *dm <= 0.0 {
                if !flagged[m] {
                    flagged[m] = true;
                    flags.push((m, FactorError::NotSpd { index: j }));
                }
                *dm = 1.0;
            }
        }
        {
            let ljj = l.plane_mut(j, j);
            for m in 0..n {
                ljj[m] = det::sqrt(d[m]);
            }
        }
        for i in (j + 1)..k {
            // L[i][j] = (A[i][j] − Σ_{p<j} L[i][p]·L[j][p]) / L[j][j]
            let mut s: Vec<f64> = a.plane(i, j).to_vec();
            for p in 0..j {
                let (lip, ljp) = (l.plane(i, p), l.plane(j, p));
                for m in 0..n {
                    s[m] = lip[m].mul_add(-ljp[m], s[m]);
                }
            }
            let ljj: Vec<f64> = l.plane(j, j).to_vec();
            let lij = l.plane_mut(i, j);
            for m in 0..n {
                lij[m] = s[m] / ljj[m];
            }
        }
    }
    (l, flags)
}

/// Batched forward substitution: solve L·y = b in place (b becomes y).
/// Lower-triangular, non-unit diagonal.
///
/// # Panics
/// On shape or batch-length mismatch.
pub fn batch_solve_lower(l: &BatchMat, b: &mut BatchVec) {
    let (k, n) = (l.k, l.n);
    assert!(b.k == k && b.n == n, "batch_solve_lower shape mismatch");
    for i in 0..k {
        let mut s: Vec<f64> = b.plane(i).to_vec();
        for j in 0..i {
            let lij = l.plane(i, j);
            let bj: Vec<f64> = b.plane(j).to_vec();
            for m in 0..n {
                s[m] = lij[m].mul_add(-bj[m], s[m]);
            }
        }
        let lii = l.plane(i, i);
        let bi = b.plane_mut(i);
        for m in 0..n {
            bi[m] = s[m] / lii[m];
        }
    }
}

/// Batched back substitution: solve U·x = b in place where U is the
/// TRANSPOSE of the given lower factor (u_ij = l_ji), i.e. the
/// Cholesky second half; or an upper factor stored explicitly when
/// `transposed_lower` is false.
pub fn batch_solve_upper(u: &BatchMat, b: &mut BatchVec, transposed_lower: bool) {
    let (k, n) = (u.k, u.n);
    assert!(b.k == k && b.n == n, "batch_solve_upper shape mismatch");
    for i in (0..k).rev() {
        let mut s: Vec<f64> = b.plane(i).to_vec();
        for j in (i + 1)..k {
            let uij = if transposed_lower {
                u.plane(j, i)
            } else {
                u.plane(i, j)
            };
            let bj: Vec<f64> = b.plane(j).to_vec();
            for m in 0..n {
                s[m] = uij[m].mul_add(-bj[m], s[m]);
            }
        }
        let uii = u.plane(i, i);
        let bi = b.plane_mut(i);
        for m in 0..n {
            bi[m] = s[m] / uii[m];
        }
    }
}

/// Convenience: solve A·x = b per matrix from a Cholesky factor
/// (forward then transposed-back substitution), in place.
pub fn batch_cholesky_solve(l: &BatchMat, b: &mut BatchVec) {
    batch_solve_lower(l, b);
    batch_solve_upper(l, b, true);
}

// ---------------------------------------------------------------------- LU

/// Batched LU with partial pivoting. Pivot selection is per-matrix
/// data-dependent (strictly-greater comparison ⇒ the LOWEST row index
/// among maximal |pivot| — `factor::lu`'s deterministic tie-break), so
/// this kernel iterates matrices in the scalar dimension; correctness
/// and determinism first, lane-width vectorization is the recorded
/// perf lane. Exactly-zero pivot columns are flagged `Singular` and
/// that matrix continues with a unit pivot (output unspecified).
pub struct BatchLu {
    /// Compact LU factors (unit-lower below, U on and above the
    /// diagonal), per matrix.
    pub lu: BatchMat,
    /// Row permutation per matrix: `perm[step * batch + m]`.
    pub perm: Vec<u32>,
    /// Per-matrix failures (empty = every matrix factored).
    pub flags: Vec<(usize, FactorError)>,
}

/// Factor every matrix in the batch: PA = LU.
///
/// Loop shape (bead 9ekv): step-OUTER. Pivot selection and the row
/// swap stay per-matrix scalar (data-dependent choices), but the
/// trailing update runs plane-wise with SIMD lanes across the batch.
/// Lanes are independent matrices and every per-element operation and
/// its order are unchanged from the matrix-at-a-time formulation, so
/// outputs are BIT-IDENTICAL (membership invariance + golden, tested).
#[must_use]
pub fn batch_lu(a: &BatchMat) -> BatchLu {
    let (k, n) = (a.k, a.n);
    let mut lu = a.clone();
    let mut perm = vec![0u32; k * n];
    let mut flags = Vec::new();
    let mut flagged = vec![false; n];
    // pstate[m·k + row]: matrix m's permutation, swapped as steps pivot.
    let kk = u32::try_from(k).expect("k fits u32");
    let mut pstate: Vec<u32> = (0..n).flat_map(|_| 0..kk).collect();
    let mut fvec = vec![0.0f64; n];
    for step in 0..k {
        // Pivot: lowest row index among maximal |value| (strict >),
        // then swap + exactly-zero flagging — per matrix. (`m` indexes
        // three parallel per-matrix structures, not one slice.)
        #[allow(clippy::needless_range_loop)]
        for m in 0..n {
            let mut best = step;
            let mut best_abs = lu.get(m, step, step).abs();
            for r in (step + 1)..k {
                let v = lu.get(m, r, step).abs();
                if v > best_abs {
                    best_abs = v;
                    best = r;
                }
            }
            if best != step {
                for c in 0..k {
                    let hi = lu.get(m, best, c);
                    let lo = lu.get(m, step, c);
                    lu.set(m, best, c, lo);
                    lu.set(m, step, c, hi);
                }
                pstate.swap(m * k + step, m * k + best);
            }
            if lu.get(m, step, step) == 0.0 {
                if !flagged[m] {
                    flagged[m] = true;
                    flags.push((m, FactorError::Singular { index: step }));
                }
                lu.set(m, step, step, 1.0);
            }
        }
        // Trailing update, lane-vectorized across the batch.
        for r in (step + 1)..k {
            {
                let (fr, piv) = lu.planes_mut2((r, step), (step, step));
                for ((fv, f), &p) in fvec.iter_mut().zip(fr.iter_mut()).zip(&*piv) {
                    *f /= p;
                    *fv = *f;
                }
            }
            for c in (step + 1)..k {
                let (rc, sc) = lu.planes_mut2((r, c), (step, c));
                for ((v, &sv), &f) in rc.iter_mut().zip(&*sc).zip(&fvec) {
                    *v = f.mul_add(-sv, *v);
                }
            }
        }
    }
    for m in 0..n {
        for step in 0..k {
            perm[step * n + m] = pstate[m * k + step];
        }
    }
    BatchLu { lu, perm, flags }
}

impl BatchLu {
    /// Solve A·x = b per matrix, in place (applies P, then L, then U).
    ///
    /// # Panics
    /// On shape or batch-length mismatch.
    pub fn solve(&self, b: &mut BatchVec) {
        let (k, n) = (self.lu.k, self.lu.n);
        assert!(b.k == k && b.n == n, "BatchLu::solve shape mismatch");
        // Apply the permutation per matrix.
        for m in 0..n {
            let x: Vec<f64> = (0..k)
                .map(|step| b.get(m, self.perm[step * n + m] as usize))
                .collect();
            for (i, v) in x.into_iter().enumerate() {
                b.set(m, i, v);
            }
        }
        // Forward (unit lower), fixed order.
        for i in 0..k {
            let mut s: Vec<f64> = b.plane(i).to_vec();
            for j in 0..i {
                let lij = self.lu.plane(i, j);
                let bj: Vec<f64> = b.plane(j).to_vec();
                for m in 0..n {
                    s[m] = lij[m].mul_add(-bj[m], s[m]);
                }
            }
            b.plane_mut(i).copy_from_slice(&s);
        }
        // Backward (U with diagonal).
        batch_solve_upper(&self.lu, b, false);
    }
}

// -------------------------------------------------------------------- eigen

/// Batched symmetric-3×3 eigenvalues, closed form (deterministic
/// trigonometric method through fs-math strict kernels), ascending.
/// The stress-principal-direction hot path.
///
/// # Panics
/// If `k != 3`.
#[must_use]
pub fn batch_eigh3_values(a: &BatchMat) -> BatchVec {
    assert_eq!(a.k, 3, "batch_eigh3_values requires k = 3");
    let n = a.n;
    let mut out = BatchVec::zeros(3, n);
    for m in 0..n {
        let g = gather3(a, m);
        let vals = eigh3_closed(&g);
        for (i, v) in vals.iter().enumerate() {
            out.set(m, i, *v);
        }
    }
    out
}

fn eigh3_closed(g: &[f64; 9]) -> [f64; 3] {
    let q = (g[0] + g[4] + g[8]) / 3.0;
    let off2 = 2.0 * g[1].mul_add(g[1], g[2].mul_add(g[2], g[5] * g[5]));
    let d0 = g[0] - q;
    let d1 = g[4] - q;
    let d2 = g[8] - q;
    let p2 = d0.mul_add(d0, d1.mul_add(d1, d2.mul_add(d2, off2)));
    if p2 == 0.0 {
        // Already diagonal-with-equal-entries: λ = q (triple).
        return [q, q, q];
    }
    let p = det::sqrt(p2 / 6.0);
    // B = (A − qI)/p; r = det(B)/2 ∈ [−1, 1] up to roundoff.
    let b = [
        (g[0] - q) / p,
        g[1] / p,
        g[2] / p,
        g[3] / p,
        (g[4] - q) / p,
        g[5] / p,
        g[6] / p,
        g[7] / p,
        (g[8] - q) / p,
    ];
    let r = (det3(&b) / 2.0).clamp(-1.0, 1.0);
    // φ = acos(r)/3 via the strict atan2/sqrt composition.
    let phi = det::atan2(det::sqrt((1.0 - r) * (1.0 + r)), r) / 3.0;
    // Roots q + 2p·cos(φ + 2πk/3), φ ∈ [0, π/3]: k = 0 is the maximum
    // (cos ∈ [½, 1]), k = 1 the minimum (cos ∈ [−1, −½]); the middle
    // one comes from the trace identity (no third cos evaluation).
    let two_pi_3 = 2.0 * std::f64::consts::FRAC_PI_3;
    let lmax = q + 2.0 * p * det::cos(phi);
    let lmin = q + 2.0 * p * det::cos(phi + two_pi_3);
    let lmid = 3.0f64.mul_add(q, -(lmax + lmin));
    [lmin, lmid, lmax]
}

/// Batched symmetric eigendecomposition for general small k (covers
/// the 6×6 metric-tensor path): per-matrix Jacobi through
/// `crate::eigen::jacobi_eigh` (exactly orthogonal vectors,
/// deterministic sweep order). Returns (values ascending, vectors) —
/// vectors[m] column j is the eigenvector of value j, stored as a
/// BatchMat.
///
/// # Panics
/// On shape mismatch between the two outputs and `a`.
#[must_use]
pub fn batch_eigh(a: &BatchMat) -> (BatchVec, BatchMat) {
    let (k, n) = (a.k, a.n);
    let mut vals = BatchVec::zeros(k, n);
    let mut vecs = BatchMat::zeros(k, n);
    for m in 0..n {
        let g = a.gather(m);
        let (w, v) = crate::eigen::jacobi_eigh(&g, k);
        for (i, wi) in w.iter().enumerate() {
            vals.set(m, i, *wi);
        }
        vecs.scatter(m, &v);
    }
    (vals, vecs)
}

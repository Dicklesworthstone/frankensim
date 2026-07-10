//! Batched small-dense f32 and MIXED-precision GEMM (bead 9ekv, scope e)
//! — the intended substrate for a future LBM moment path: distributions
//! stored f32 (half the plane traffic of f64), moments accumulated in
//! f64 (widening is EXACT for every f32), narrowed ONCE at the store.
//! There is no production `fs-lbm` consumer yet. Same entry-plane SoA
//! layout and per-element determinism doctrine as [`crate::batched`]:
//!
//! - one fused operation order per element (l-ascending `mul_add`),
//!   identical for every lane — batch membership can never change bits;
//! - β = 0 overwrites (BLAS convention; NaN in C is overwritten);
//! - the mixed kernel's contract is "f64 chain + exactly one f32
//!   rounding per output element", which is bit-deterministic on every
//!   conforming target.
//!
//! v1 is the CORRECTNESS layer: plain plane sweeps, no chunking or
//! capsule dispatch — the perf treatment (MBLK partitioning, packed
//! tiles) follows the f64 path's 9ekv lane once these semantics are
//! consumed (no performance claim yet; see the crate contract).

use fs_soa::FieldBuf;

/// Semantic version of the batched-f32/mixed bit contract
/// (golden-couplings surface `fs-la:batched-f32-bits`): the fused
/// per-element chain order, the widen-exact/narrow-once mixed rule,
/// β = 0 overwrite, and α = 0 no-read convention. Changing any of them bumps this and
/// deliberately re-freezes the dependents in golden-couplings.json.
pub const BATCHED_F32_BIT_SEMANTICS_VERSION: u32 = 2;

/// Plane stride granularity: 32 f32 = 128 bytes, so every plane start
/// stays 128-byte aligned given the buffer base is.
const STRIDE_QUANTUM: usize = 32;

fn stride_for(n: usize) -> usize {
    let blocks = n / STRIDE_QUANTUM + usize::from(!n.is_multiple_of(STRIDE_QUANTUM));
    let s = blocks
        .checked_mul(STRIDE_QUANTUM)
        .expect("batch stride overflow");
    // Break power-of-two plane BYTE spacing (the f64 path measured 2×
    // set-aliasing slowdowns; 8192 f32 = the same 32 KB byte spacing
    // where its rule engages). Pure layout — bit-neutral.
    if s >= 8192 && s.is_power_of_two() {
        s.checked_add(STRIDE_QUANTUM)
            .expect("batch stride padding overflow")
    } else {
        s
    }
}

fn make_buf(planes: usize, stride: usize) -> FieldBuf<f32> {
    let total = planes
        .checked_mul(stride)
        .expect("batch allocation size overflow");
    let mut data = FieldBuf::with_capacity(total);
    for _ in 0..total {
        data.push(0.0f32);
    }
    data
}

/// A batch of `n` dense k×k f32 matrices in entry-plane SoA layout:
/// plane p = i·k + j holds entry (i, j) of every matrix, 128-byte
/// aligned.
#[derive(Debug, Clone)]
pub struct BatchMatF32 {
    data: FieldBuf<f32>,
    k: usize,
    n: usize,
    stride: usize,
}

impl BatchMatF32 {
    /// Batch of `n` zero k×k matrices.
    ///
    /// # Panics
    /// If either dimension is zero or the padded allocation shape overflows.
    #[must_use]
    pub fn zeros(k: usize, n: usize) -> BatchMatF32 {
        assert!(k > 0 && n > 0, "empty batch shapes are programmer errors");
        let stride = stride_for(n);
        let planes = k.checked_mul(k).expect("batch matrix shape overflow");
        BatchMatF32 {
            data: make_buf(planes, stride),
            k,
            n,
            stride,
        }
    }

    /// Batch filled from `f(m, i, j)` (matrix index, row, column).
    #[must_use]
    pub fn from_fn<F: FnMut(usize, usize, usize) -> f32>(
        k: usize,
        n: usize,
        mut f: F,
    ) -> BatchMatF32 {
        let mut out = BatchMatF32::zeros(k, n);
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
    pub fn get(&self, m: usize, i: usize, j: usize) -> f32 {
        self.plane(i, j)[m]
    }

    /// Set entry (i, j) of matrix `m`.
    pub fn set(&mut self, m: usize, i: usize, j: usize, v: f32) {
        self.plane_mut(i, j)[m] = v;
    }

    /// The (i, j) entry plane across the whole batch (len = batch).
    #[must_use]
    pub fn plane(&self, i: usize, j: usize) -> &[f32] {
        assert!(i < self.k && j < self.k, "plane index out of bounds");
        let p = i * self.k + j;
        &self.data.as_slice()[p * self.stride..p * self.stride + self.n]
    }

    /// Mutable (i, j) entry plane across the whole batch.
    pub fn plane_mut(&mut self, i: usize, j: usize) -> &mut [f32] {
        assert!(i < self.k && j < self.k, "plane index out of bounds");
        let p = i * self.k + j;
        let (stride, n) = (self.stride, self.n);
        &mut self.data.as_mut_slice()[p * stride..p * stride + n]
    }
}

/// Pure-f32 batched GEMM: C ← α·A·B + β·C per matrix, one fused f32
/// `mul_add` chain per output element (l ascending), β = 0 overwrites.
/// Bit-deterministic and batch-membership invariant by construction.
///
/// # Panics
/// On dimension or batch-length mismatches (programmer errors).
pub fn batch_gemm_f32(
    alpha: f32,
    a: &BatchMatF32,
    b: &BatchMatF32,
    beta: f32,
    c: &mut BatchMatF32,
) {
    let k = a.k;
    assert!(b.k == k && c.k == k, "batch_gemm_f32 dimension mismatch");
    assert!(
        a.n == b.n && a.n == c.n,
        "batch_gemm_f32 batch-length mismatch"
    );
    if alpha == 0.0 {
        scale_batch_f32(beta, c);
        return;
    }
    let n = a.n;
    let mut acc = vec![0.0f32; n];
    for i in 0..k {
        for j in 0..k {
            acc.fill(0.0);
            for l in 0..k {
                let ap = a.plane(i, l);
                let bp = b.plane(l, j);
                for ((s, &am), &bm) in acc.iter_mut().zip(ap).zip(bp) {
                    *s = am.mul_add(bm, *s);
                }
            }
            let cp = c.plane_mut(i, j);
            if beta == 0.0 {
                for (cm, &s) in cp.iter_mut().zip(&acc) {
                    *cm = alpha * s;
                }
            } else {
                for (cm, &s) in cp.iter_mut().zip(&acc) {
                    *cm = alpha.mul_add(s, beta * *cm);
                }
            }
        }
    }
}

/// MIXED batched GEMM for the LBM moment path: operands and result
/// stored f32, the whole per-element computation carried in f64 —
/// widening every f32 is EXACT, the l-ascending chain is fused f64
/// `mul_add`, α/β apply in f64 (β = 0 overwrites; C's old value widens
/// exactly), and the result narrows to f32 exactly ONCE. Contract:
/// each output element experiences a single f32 rounding.
///
/// # Panics
/// On dimension or batch-length mismatches (programmer errors).
#[allow(clippy::cast_possible_truncation)] // the narrow IS the contract
pub fn batch_gemm_mixed(
    alpha: f64,
    a: &BatchMatF32,
    b: &BatchMatF32,
    beta: f64,
    c: &mut BatchMatF32,
) {
    let k = a.k;
    assert!(b.k == k && c.k == k, "batch_gemm_mixed dimension mismatch");
    assert!(
        a.n == b.n && a.n == c.n,
        "batch_gemm_mixed batch-length mismatch"
    );
    if alpha == 0.0 {
        scale_batch_mixed(beta, c);
        return;
    }
    let n = a.n;
    let mut acc = vec![0.0f64; n];
    for i in 0..k {
        for j in 0..k {
            acc.fill(0.0);
            for l in 0..k {
                let ap = a.plane(i, l);
                let bp = b.plane(l, j);
                for ((s, &am), &bm) in acc.iter_mut().zip(ap).zip(bp) {
                    *s = f64::from(am).mul_add(f64::from(bm), *s);
                }
            }
            let cp = c.plane_mut(i, j);
            if beta == 0.0 {
                for (cm, &s) in cp.iter_mut().zip(&acc) {
                    *cm = (alpha * s) as f32;
                }
            } else {
                for (cm, &s) in cp.iter_mut().zip(&acc) {
                    *cm = alpha.mul_add(s, beta * f64::from(*cm)) as f32;
                }
            }
        }
    }
}

fn scale_batch_f32(beta: f32, c: &mut BatchMatF32) {
    for i in 0..c.k {
        for j in 0..c.k {
            let cp = c.plane_mut(i, j);
            if beta == 0.0 {
                cp.fill(0.0);
            } else if beta != 1.0 {
                for value in cp {
                    *value *= beta;
                }
            }
        }
    }
}

#[allow(clippy::cast_possible_truncation)] // one narrowing is the mixed contract
fn scale_batch_mixed(beta: f64, c: &mut BatchMatF32) {
    for i in 0..c.k {
        for j in 0..c.k {
            let cp = c.plane_mut(i, j);
            if beta == 0.0 {
                cp.fill(0.0);
            } else if beta != 1.0 {
                for value in cp {
                    *value = (beta * f64::from(*value)) as f32;
                }
            }
        }
    }
}

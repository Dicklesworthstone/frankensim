//! BLIS-style GEMM (plan §6.1): C = α·A·B + β·C, row-major, with cache
//! blocking and panel packing. v1 is the CORRECTNESS + ARCHITECTURE layer
//! in safe Rust: the register-tiled microkernel accumulates in local
//! arrays with fused `mul_add` (auto-vectorizes respectably); the
//! arch-specific fs-simd capsule microkernels and the autotuned blocking
//! sweep are the recorded perf follow-up (gated on the autotuner bead).
//!
//! DETERMINISM CONTRACT: accumulation order is k-ascending within KC
//! chunks, with per-chunk register partials folded into C in chunk order.
//! Hence KC is PART of the bit contract (changing it legitimately changes
//! bits → golden bump with justification); MC/NC are bit-neutral (pure
//! m/n tiling — the fact the future parallel lane relies on). Everything
//! is fixed-order +/×/mul_add: cross-ISA bit-deterministic by
//! construction, golden-hashed in tests.

/// Micro-tile rows (A panel height). Pre-autotuner default.
const MR: usize = 8;
/// Micro-tile cols (B panel width). Pre-autotuner default.
const NR: usize = 4;
/// K blocking — PART OF THE BIT CONTRACT (see module docs).
const KC: usize = 256;
/// M blocking (bit-neutral).
const MC: usize = 128;
/// N blocking (bit-neutral).
const NC: usize = 512;

/// f64 GEMM: `c[m×n] = alpha · a[m×k] · b[k×n] + beta · c`, row-major
/// contiguous slices. β = 0 OVERWRITES c (existing NaN/garbage in c is
/// ignored — the BLAS convention callers expect for uninitialized output).
///
/// # Panics
/// Structured panics on slice-length mismatches.
#[allow(clippy::too_many_arguments)] // BLAS-shape signature: m,n,k,alpha,a,b,beta,c
pub fn gemm_f64(
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
) {
    assert_eq!(a.len(), m * k, "a must be m*k = {}", m * k);
    assert_eq!(b.len(), k * n, "b must be k*n = {}", k * n);
    assert_eq!(c.len(), m * n, "c must be m*n = {}", m * n);
    // β pass first (once, before any KC chunk): scale or overwrite.
    scale_c(c, beta);
    if m == 0 || n == 0 || alpha == 0.0 {
        return;
    }
    if k == 0 {
        return; // C = beta*C only (already applied).
    }
    let mut a_pack = vec![0.0f64; MC * KC];
    let mut b_pack = vec![0.0f64; KC * NC];
    // Loop nest (BLIS order): NC → KC → MC → NR → MR → K.
    let mut jc = 0;
    while jc < n {
        let nc = NC.min(n - jc);
        let mut pc = 0;
        while pc < k {
            let kc = KC.min(k - pc);
            pack_b(&mut b_pack, b, n, pc, jc, kc, nc);
            let mut ic = 0;
            while ic < m {
                let mc = MC.min(m - ic);
                pack_a(&mut a_pack, a, k, ic, pc, mc, kc);
                macro_kernel(&a_pack, &b_pack, c, m, n, ic, jc, mc, nc, kc, alpha);
                ic += MC;
            }
            pc += KC;
        }
        jc += NC;
    }
    let _ = m; // (m used above; silences pedantic when MC >= m)
}

/// β application with BLAS overwrite semantics for β = 0.
fn scale_c(c: &mut [f64], beta: f64) {
    if beta == 0.0 {
        c.fill(0.0);
    } else if beta.to_bits() != 1.0f64.to_bits() {
        for v in c.iter_mut() {
            *v *= beta;
        }
    }
}

/// Pack an mc×kc block of A (row-major, ld = k) into MR-row micro-panels:
/// panel p holds rows [p·MR, p·MR+MR) column-major-within-panel
/// (k-index fastest across the MR lanes). Short tail rows are zero-padded —
/// zero lanes contribute exact +0.0 products which never reach C (tail
/// handling masks them on write-back).
fn pack_a(dst: &mut [f64], a: &[f64], lda: usize, ic: usize, pc: usize, mc: usize, kc: usize) {
    let mut w = 0;
    let mut p = 0;
    while p < mc {
        let rows = MR.min(mc - p);
        for kk in 0..kc {
            for r in 0..MR {
                dst[w] = if r < rows {
                    a[(ic + p + r) * lda + pc + kk]
                } else {
                    0.0
                };
                w += 1;
            }
        }
        p += MR;
    }
}

/// Pack a kc×nc block of B (row-major, ld = n) into NR-column micro-panels
/// (k-index outer, NR lanes inner), zero-padded tails.
fn pack_b(dst: &mut [f64], b: &[f64], ldb: usize, pc: usize, jc: usize, kc: usize, nc: usize) {
    let mut w = 0;
    let mut q = 0;
    while q < nc {
        let cols = NR.min(nc - q);
        for kk in 0..kc {
            for s in 0..NR {
                dst[w] = if s < cols {
                    b[(pc + kk) * ldb + jc + q + s]
                } else {
                    0.0
                };
                w += 1;
            }
        }
        q += NR;
    }
}

/// The macro kernel: sweep micro-tiles of the packed panels.
#[allow(clippy::too_many_arguments)]
fn macro_kernel(
    a_pack: &[f64],
    b_pack: &[f64],
    c: &mut [f64],
    _m: usize,
    n: usize,
    ic: usize,
    jc: usize,
    mc: usize,
    nc: usize,
    kc: usize,
    alpha: f64,
) {
    let mut p = 0;
    while p < mc {
        let rows = MR.min(mc - p);
        let a_panel = &a_pack[(p / MR) * MR * kc..][..MR * kc];
        let mut q = 0;
        while q < nc {
            let cols = NR.min(nc - q);
            let b_panel = &b_pack[(q / NR) * NR * kc..][..NR * kc];
            // Register-tiled microkernel: MR×NR accumulators, k ascending.
            let mut acc = [[0.0f64; NR]; MR];
            for kk in 0..kc {
                let av = &a_panel[kk * MR..kk * MR + MR];
                let bv = &b_panel[kk * NR..kk * NR + NR];
                for (accr, &ar) in acc.iter_mut().zip(av) {
                    for (slot, &bs) in accr.iter_mut().zip(bv) {
                        *slot = ar.mul_add(bs, *slot);
                    }
                }
            }
            // Write-back with α, masking padded tail lanes.
            for (r, accr) in acc.iter().enumerate().take(rows) {
                let crow = (ic + p + r) * n + jc + q;
                for (s, &av) in accr.iter().enumerate().take(cols) {
                    c[crow + s] = alpha.mul_add(av, c[crow + s]);
                }
            }
            q += NR;
        }
        p += MR;
    }
}

/// f32 GEMM, same structure and contract (KC shared).
#[allow(clippy::too_many_arguments)] // BLAS-shape signature: m,n,k,alpha,a,b,beta,c
pub fn gemm_f32(
    m: usize,
    n: usize,
    k: usize,
    alpha: f32,
    a: &[f32],
    b: &[f32],
    beta: f32,
    c: &mut [f32],
) {
    assert_eq!(a.len(), m * k, "a must be m*k = {}", m * k);
    assert_eq!(b.len(), k * n, "b must be k*n = {}", k * n);
    assert_eq!(c.len(), m * n, "c must be m*n = {}", m * n);
    if beta == 0.0 {
        c.fill(0.0);
    } else if beta.to_bits() != 1.0f32.to_bits() {
        for v in c.iter_mut() {
            *v *= beta;
        }
    }
    if m == 0 || n == 0 || k == 0 || alpha == 0.0 {
        return;
    }
    // v1 f32 path: same loop order, unpacked (packing joins the perf bead;
    // correctness and determinism contracts are identical).
    let mut pc = 0;
    while pc < k {
        let kc = KC.min(k - pc);
        for i in 0..m {
            for j in 0..n {
                let mut acc = 0.0f32;
                for kk in 0..kc {
                    acc = a[i * k + pc + kk].mul_add(b[(pc + kk) * n + j], acc);
                }
                c[i * n + j] = alpha.mul_add(acc, c[i * n + j]);
            }
        }
        pc += KC;
    }
}

/// Mixed-precision GEMM: f32 STORAGE, f64 ACCUMULATION — the
/// bandwidth-vs-accuracy compromise used throughout the plan (§6.1). Each
/// f32 element is widened exactly (f32→f64 is exact); all arithmetic is
/// f64 mul_add in the same k-ascending order.
#[allow(clippy::too_many_arguments)] // BLAS-shape signature: m,n,k,alpha,a,b,beta,c
pub fn gemm_mixed(
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f32],
    b: &[f32],
    beta: f64,
    c: &mut [f64],
) {
    assert_eq!(a.len(), m * k, "a must be m*k = {}", m * k);
    assert_eq!(b.len(), k * n, "b must be k*n = {}", k * n);
    assert_eq!(c.len(), m * n, "c must be m*n = {}", m * n);
    scale_c(c, beta);
    if m == 0 || n == 0 || k == 0 || alpha == 0.0 {
        return;
    }
    let mut pc = 0;
    while pc < k {
        let kc = KC.min(k - pc);
        for i in 0..m {
            for j in 0..n {
                let mut acc = 0.0f64;
                for kk in 0..kc {
                    let av = f64::from(a[i * k + pc + kk]);
                    let bv = f64::from(b[(pc + kk) * n + j]);
                    acc = av.mul_add(bv, acc);
                }
                c[i * n + j] = alpha.mul_add(acc, c[i * n + j]);
            }
        }
        pc += KC;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    }

    /// The oracle: naive triple loop with the SAME KC chunking and fused
    /// arithmetic — bitwise-comparable; plus a plain tolerance oracle.
    #[allow(clippy::too_many_arguments)] // mirrors the BLAS-shape signature
    fn naive_chunked(
        m: usize,
        n: usize,
        k: usize,
        alpha: f64,
        a: &[f64],
        b: &[f64],
        beta: f64,
        c0: &[f64],
    ) -> Vec<f64> {
        let mut c: Vec<f64> = if beta == 0.0 {
            vec![0.0; m * n]
        } else {
            c0.iter()
                .map(|&v| {
                    if beta.to_bits() == 1.0f64.to_bits() {
                        v
                    } else {
                        v * beta
                    }
                })
                .collect()
        };
        let mut pc = 0;
        while pc < k {
            let kc = KC.min(k - pc);
            for i in 0..m {
                for j in 0..n {
                    let mut acc = 0.0f64;
                    for kk in 0..kc {
                        acc = a[i * k + pc + kk].mul_add(b[(pc + kk) * n + j], acc);
                    }
                    c[i * n + j] = alpha.mul_add(acc, c[i * n + j]);
                }
            }
            pc += KC;
        }
        c
    }

    fn rand_mat(rows: usize, cols: usize, seed: u64) -> Vec<f64> {
        let mut s = seed;
        (0..rows * cols).map(|_| lcg(&mut s)).collect()
    }

    #[test]
    fn matches_oracle_bitwise_across_shape_sweep() {
        // The packed/blocked path must be BIT-IDENTICAL to the same-order
        // naive path: packing must not change arithmetic, only layout.
        let shapes = [
            (1usize, 1usize, 1usize),
            (1, 7, 3),
            (5, 1, 9),
            (8, 4, 256),   // exactly one micro-tile, one KC chunk
            (9, 5, 257),   // tails in every dimension
            (33, 17, 300), // KC chunking engaged
            (64, 64, 64),
            (3, 200, 2), // wide
            (200, 3, 2), // tall-skinny
        ];
        for (idx, &(m, n, k)) in shapes.iter().enumerate() {
            let a = rand_mat(m, k, 0xA + idx as u64);
            let b = rand_mat(k, n, 0xB + idx as u64);
            let c0 = rand_mat(m, n, 0xC + idx as u64);
            for (alpha, beta) in [(1.0, 0.0), (2.5, 1.0), (-0.75, 0.5)] {
                let mut c = c0.clone();
                gemm_f64(m, n, k, alpha, &a, &b, beta, &mut c);
                let want = naive_chunked(m, n, k, alpha, &a, &b, beta, &c0);
                for (i, (&got, &w)) in c.iter().zip(&want).enumerate() {
                    assert_eq!(
                        got.to_bits(),
                        w.to_bits(),
                        "({m}x{n}x{k}) alpha={alpha} beta={beta} at {i}: {got} vs {w}"
                    );
                }
            }
        }
        println!(
            "{{\"suite\":\"fs-la\",\"case\":\"gemm-oracle\",\"verdict\":\"pass\",\"detail\":\"9 shapes x 3 (alpha,beta) bitwise vs same-order oracle\"}}"
        );
    }

    #[test]
    fn degenerate_and_beta_semantics() {
        // k = 0: C = beta*C, and beta = 0 must OVERWRITE garbage/NaN.
        let mut c = vec![f64::NAN, 3.0, -2.0, 1.0];
        gemm_f64(2, 2, 0, 1.0, &[], &[], 0.0, &mut c);
        assert!(
            c.iter().all(|&v| v == 0.0),
            "beta=0 must overwrite NaN: {c:?}"
        );
        let mut c2 = vec![1.0, 2.0, 3.0, 4.0];
        gemm_f64(2, 2, 0, 1.0, &[], &[], 2.0, &mut c2);
        assert_eq!(c2, vec![2.0, 4.0, 6.0, 8.0]);
        // m or n zero: no-op, no panic.
        let mut empty: Vec<f64> = vec![];
        gemm_f64(0, 5, 3, 1.0, &[], &rand_mat(3, 5, 1), 0.0, &mut empty);
        // alpha = 0 leaves beta*C.
        let a = rand_mat(2, 3, 2);
        let b = rand_mat(3, 2, 3);
        let mut c3 = vec![1.0; 4];
        gemm_f64(2, 2, 3, 0.0, &a, &b, 1.0, &mut c3);
        assert_eq!(c3, vec![1.0; 4]);
    }

    #[test]
    fn transpose_identity_and_submatrix_consistency() {
        let (m, n, k) = (24usize, 18usize, 40usize);
        let a = rand_mat(m, k, 7);
        let b = rand_mat(k, n, 8);
        // (A·B)ᵀ == Bᵀ·Aᵀ within tight tolerance (orders differ → not bitwise).
        let mut ab = vec![0.0; m * n];
        gemm_f64(m, n, k, 1.0, &a, &b, 0.0, &mut ab);
        let at: Vec<f64> = (0..k * m).map(|i| a[(i % m) * k + i / m]).collect();
        let bt: Vec<f64> = (0..n * k).map(|i| b[(i % k) * n + i / k]).collect();
        let mut btat = vec![0.0; n * m];
        gemm_f64(n, m, k, 1.0, &bt, &at, 0.0, &mut btat);
        for i in 0..m {
            for j in 0..n {
                let x = ab[i * n + j];
                let y = btat[j * m + i];
                assert!(
                    (x - y).abs() <= 1e-13 * x.abs().max(1.0),
                    "transpose identity at ({i},{j}): {x} vs {y}"
                );
            }
        }
        // Submatrix consistency: the top-left quadrant of C equals the GEMM
        // of the corresponding A rows with B (exact: row tiling is
        // bit-neutral, same k order).
        let m2 = m / 2;
        let a_top = &a[..m2 * k];
        let mut c_top = vec![0.0; m2 * n];
        gemm_f64(m2, n, k, 1.0, a_top, &b, 0.0, &mut c_top);
        for i in 0..m2 * n {
            assert_eq!(
                c_top[i].to_bits(),
                ab[i].to_bits(),
                "row-tiling changed bits at {i}"
            );
        }
    }

    #[test]
    fn f32_and_mixed_paths() {
        let (m, n, k) = (17usize, 13usize, 129usize);
        let mut s = 0x32_u64;
        let af: Vec<f32> = (0..m * k).map(|_| lcg(&mut s) as f32).collect();
        let bf: Vec<f32> = (0..k * n).map(|_| lcg(&mut s) as f32).collect();
        // Mixed vs full-f64 reference on the widened inputs: mixed IS the
        // f64 computation on exactly-widened values — bitwise equal.
        let ad: Vec<f64> = af.iter().map(|&v| f64::from(v)).collect();
        let bd: Vec<f64> = bf.iter().map(|&v| f64::from(v)).collect();
        let mut c_mixed = vec![0.0f64; m * n];
        gemm_mixed(m, n, k, 1.0, &af, &bf, 0.0, &mut c_mixed);
        let c_ref = naive_chunked(m, n, k, 1.0, &ad, &bd, 0.0, &vec![0.0; m * n]);
        for i in 0..m * n {
            assert_eq!(
                c_mixed[i].to_bits(),
                c_ref[i].to_bits(),
                "mixed != widened f64 at {i}"
            );
        }
        // f32 path: tolerance vs f64 reference (f32 accumulates in f32).
        let mut c32 = vec![0.0f32; m * n];
        gemm_f32(m, n, k, 1.0, &af, &bf, 0.0, &mut c32);
        for i in 0..m * n {
            let err = (f64::from(c32[i]) - c_ref[i]).abs();
            assert!(
                err <= 1e-4 * c_ref[i].abs().max(1.0),
                "f32 path error {err} at {i}"
            );
        }
        println!(
            "{{\"suite\":\"fs-la\",\"case\":\"gemm-precisions\",\"verdict\":\"pass\",\"detail\":\"mixed == widened-f64 bitwise; f32 within 1e-4\"}}"
        );
    }

    #[test]
    fn deterministic_golden_hash() {
        let (m, n, k) = (48usize, 36usize, 300usize);
        let a = rand_mat(m, k, 0x60);
        let b = rand_mat(k, n, 0x61);
        let run = || {
            let mut c = vec![0.0; m * n];
            gemm_f64(m, n, k, 1.25, &a, &b, 0.0, &mut c);
            let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
            for v in &c {
                for byte in v.to_bits().to_le_bytes() {
                    acc ^= u64::from(byte);
                    acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
                }
            }
            acc
        };
        let h = run();
        assert_eq!(h, run(), "same inputs must give identical bits");
        println!(
            "{{\"suite\":\"fs-la\",\"case\":\"gemm-golden\",\"verdict\":\"info\",\"detail\":\"{h:#018x}\"}}"
        );
        assert_eq!(
            h, GOLDEN_HASH,
            "GEMM output bits changed: {h:#018x} vs {GOLDEN_HASH:#018x} — KC is part of the \
             bit contract; bump only with semantic justification"
        );
    }

    /// Recorded on aarch64-apple (M4 Pro); must match on x86-64 (trj).
    const GOLDEN_HASH: u64 = 0x1d7a_a3c6_b631_7ef0;
}

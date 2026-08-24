//! Cross-ISA determinism goldens for the `fs-modal` certified spectral
//! slice (bead `frankensim-music-v8-root-3ez8g.13.4`).
//!
//! `slice_window`'s inertia counting, factorizations, Lanczos iterations,
//! and deflation restarts must produce bit-identical certified mode
//! intervals on both reference ISA families in both build modes.
//!
//! Fixture discipline: the 49x49 grid Laplacian has a degenerate spectrum
//! (the (i,j)/(j,i) pairs), so the windows below are FROZEN as bit-exact
//! hex f64 endpoints — authored once offline strictly between distinct
//! consecutive analytic eigenvalues, never recomputed from libm trig at
//! test time. A cross-ISA difference in the fixture would masquerade as a
//! kernel difference; freezing the endpoints removes that failure mode.

use fs_modal::{SliceOptions, slice_window};
use fs_sparse::{Coo, Csr};

const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

fn fold(acc: u64, v: f64) -> u64 {
    v.to_bits()
        .to_le_bytes()
        .iter()
        .fold(acc, |a, &b| (a ^ u64::from(b)).wrapping_mul(FNV_PRIME))
}

fn fold_u64(acc: u64, v: u64) -> u64 {
    (acc ^ v).wrapping_mul(FNV_PRIME)
}

fn laplacian_2d(s: usize) -> Csr {
    let n = s * s;
    let mut coo = Coo::new(n, n);
    for i in 0..s {
        for j in 0..s {
            let u = i * s + j;
            coo.push(u, u, 4.0);
            if i > 0 {
                coo.push(u, u - s, -1.0);
            }
            if i + 1 < s {
                coo.push(u, u + s, -1.0);
            }
            if j > 0 {
                coo.push(u, u - 1, -1.0);
            }
            if j + 1 < s {
                coo.push(u, u + 1, -1.0);
            }
        }
    }
    coo.assemble()
}

fn identity(n: usize) -> Csr {
    let mut coo = Coo::new(n, n);
    for i in 0..n {
        coo.push(i, i, 1.0);
    }
    coo.assemble()
}

/// Window endpoints: midpoints between distinct consecutive analytic
/// eigenvalues of the 49x49 grid Laplacian, frozen as exact bit patterns
/// (`2 − 2·cos(iπ/50)` family; windows cover 9 and 8 modes respectively,
/// spanning degenerate pairs so the deflation-restart path is exercised).
const W1_LO_BITS: u64 = 0x3fa2_275e_af4f_ffd0;
const W1_HI_BITS: u64 = 0x3fb6_a376_29ce_f9f0;
const W2_LO_BITS: u64 = 0x3fc0_9320_19d9_8d08;
const W2_HI_BITS: u64 = 0x3fc5_8abe_165c_9cf4;

/// Certified counts inside each frozen window (analytic, verified at
/// authoring time against the corner-block scan used by the casebook).
const W1_EXPECTED: usize = 9;
const W2_EXPECTED: usize = 8;

fn window_digest(k: &Csr, m: &Csr, lo: f64, hi: f64, expected: usize) -> u64 {
    let rep = slice_window(k, m, (lo, hi), &SliceOptions::default()).expect("slice");
    assert_eq!(rep.modes.len(), expected, "inertia count must equal the frozen count");
    let mut acc = fold_u64(FNV_OFFSET, u64::try_from(expected).expect("count"));
    for mode in &rep.modes {
        acc = fold(acc, mode.interval.0);
        acc = fold(acc, mode.interval.1);
        acc = fold(acc, mode.residual);
    }
    acc
}

/// Verified bit-identical aarch64-apple and x86_64-linux (debug) on
/// 2026-08-23, bead frankensim-music-v8-root-3ez8g.13.4.
const GOLDEN_HASH: u64 = 0x11dc_ed94_c7f6_7115;

#[test]
fn modal_slice_digest_is_cross_isa_golden() {
    let s = 49;
    let k = laplacian_2d(s);
    let m = identity(s * s);
    let d1 = window_digest(&k, &m, f64::from_bits(W1_LO_BITS), f64::from_bits(W1_HI_BITS), W1_EXPECTED);
    let d2 = window_digest(&k, &m, f64::from_bits(W2_LO_BITS), f64::from_bits(W2_HI_BITS), W2_EXPECTED);
    let acc = fold_u64(d1, d2);

    println!(
        "{{\"suite\":\"fs-modal\",\"case\":\"cross-isa-slice\",\"arch\":\"{}\",\
         \"profile\":\"{}\",\"digest\":\"{acc:#018x}\",\"verdict\":\"golden-check\"}}",
        std::env::consts::ARCH,
        if cfg!(debug_assertions) { "debug" } else { "release" },
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "spectral slice bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — cross-ISA \
         golden event: bisect stage-wise, name the hazard, route through det:: in the \
         same commit"
    );
}

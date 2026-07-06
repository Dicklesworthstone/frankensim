//! fs-sparse conformance: the cross-format, cross-ISA battery any
//! reimplementation must pass (plan §13.3). Builds FEM-patterned and
//! adversarial matrices, runs every format's SpMV plus the pattern algebra,
//! and folds all output bits into one FNV-64 golden hash — recorded on
//! aarch64-apple and required to match on x86-64 (the same evidence
//! discipline as fs-math/fs-fft).

use fs_sparse::{Bsr, Coo, Csr, Sell, ops};

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
}

fn laplacian_2d(n: usize) -> Csr {
    let dim = n * n;
    let mut coo = Coo::new(dim, dim);
    for i in 0..n {
        for j in 0..n {
            let u = i * n + j;
            coo.push(u, u, 4.0);
            if i > 0 {
                coo.push(u, u - n, -1.0);
            }
            if i + 1 < n {
                coo.push(u, u + n, -1.0);
            }
            if j > 0 {
                coo.push(u, u - 1, -1.0);
            }
            if j + 1 < n {
                coo.push(u, u + 1, -1.0);
            }
        }
    }
    coo.assemble()
}

/// Recorded on aarch64-apple (M4 Pro); must match on x86-64 (trj) — the
/// cross-ISA determinism evidence for assembly + all SpMV kernels + SpGEMM.
const GOLDEN_HASH: u64 = 0xbcf5_52b6_c5bf_aed6;

#[test]
fn cross_format_battery_and_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for b in v.to_bits().to_le_bytes() {
            acc ^= u64::from(b);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };

    // Matrix zoo: FEM Laplacian, random rectangular, skewed (dense row +
    // empties), block-structured.
    let lap = laplacian_2d(12); // 144x144
    let mut seed = 0xC0FFEE_u64;
    let mut rand_m = Coo::new(96, 96);
    for r in 0..96 {
        for _ in 0..6 {
            let c = ((lcg(&mut seed) + 0.5) * 96.0) as usize % 96;
            rand_m.push(r, c, lcg(&mut seed));
        }
    }
    let rnd = rand_m.assemble();
    let mut skew_m = Coo::new(64, 64);
    for c in 0..64 {
        skew_m.push(20, c, 0.5 - c as f64 / 64.0);
    }
    for r in 0..64 {
        if r % 3 == 0 {
            skew_m.push(r, r, 2.0);
        }
    }
    let skew = skew_m.assemble();

    for (name, a) in [("laplacian", &lap), ("random", &rnd), ("skew", &skew)] {
        let x: Vec<f64> = (0..a.ncols()).map(|_| lcg(&mut seed)).collect();
        let mut y_csr = vec![0.0; a.nrows()];
        a.spmv(&x, &mut y_csr);

        // Every format must agree BITWISE.
        let sell = Sell::from_csr(a, 8, 32);
        let mut y_sell = vec![0.0; a.nrows()];
        sell.spmv(&x, &mut y_sell);
        for r in 0..a.nrows() {
            assert_eq!(
                y_csr[r].to_bits(),
                y_sell[r].to_bits(),
                "{name}: SELL diverged from CSR at row {r}"
            );
        }
        if a.nrows().is_multiple_of(4) && a.ncols().is_multiple_of(4) {
            let bsr = Bsr::from_csr(a, 4, 4);
            let mut y_bsr = vec![0.0; a.nrows()];
            bsr.spmv(&x, &mut y_bsr);
            for r in 0..a.nrows() {
                assert_eq!(
                    y_csr[r].to_bits(),
                    y_bsr[r].to_bits(),
                    "{name}: BSR diverged from CSR at row {r}"
                );
            }
        }
        for &v in &y_csr {
            feed(v);
        }

        // Pattern algebra folded in: transpose SpMV, symmetrized SpMV, A·Aᵀ.
        let at = ops::transpose(a);
        let mut y_t = vec![0.0; at.nrows()];
        let xt: Vec<f64> = (0..at.ncols()).map(|_| lcg(&mut seed)).collect();
        at.spmv(&xt, &mut y_t);
        for &v in &y_t {
            feed(v);
        }
        if a.nrows() == a.ncols() {
            let s = ops::symmetrize(a);
            let mut y_s = vec![0.0; s.nrows()];
            s.spmv(&x, &mut y_s);
            for &v in &y_s {
                feed(v);
            }
        }
        let aat = ops::spgemm(a, &at);
        let mut y_g = vec![0.0; aat.nrows()];
        let xg: Vec<f64> = (0..aat.ncols()).map(|_| lcg(&mut seed)).collect();
        aat.spmv(&xg, &mut y_g);
        for &v in &y_g {
            feed(v);
        }
    }

    println!(
        "{{\"suite\":\"fs-sparse\",\"case\":\"golden-hash\",\"verdict\":\"info\",\"detail\":\"{acc:#018x}\"}}"
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "sparse kernel output bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only \
         with semantic justification (golden-evidence policy)"
    );
}

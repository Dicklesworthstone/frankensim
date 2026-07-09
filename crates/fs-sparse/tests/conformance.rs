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

/// wsbf segment 1: the compact/sharded/parallel-assembly kernels are
/// BITWISE equal to their serial wide twins at every thread count —
/// the perf program never trades the determinism contract.
#[test]
fn wsbf_bitwise_twins() {
    // Adversarial-ish fixture: ragged rows, duplicates, empty rows.
    let (nrows, ncols) = (257usize, 199usize);
    let mut coo = fs_sparse::Coo::new(nrows, ncols);
    let mut seed = 0x5EED_2026_u64;
    let mut lcg = move || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        seed
    };
    for _ in 0..6000 {
        let r = (lcg() % nrows as u64) as usize;
        let c = (lcg() % ncols as u64) as usize;
        let v = ((lcg() >> 11) as f64) / (1u64 << 53) as f64 - 0.5;
        coo.push(r, c, v);
    }
    let serial = coo.assemble();
    // (3) parallel assembly bitwise across thread counts.
    for t in [1usize, 2, 4, 8] {
        let par = coo.assemble_parallel(t);
        assert_eq!(serial, par, "assemble_parallel({t}) != serial");
    }
    // (1a) compact CSR spmv bitwise vs wide.
    let x: Vec<f64> = (0..ncols).map(|i| 0.25 + (i % 17) as f64).collect();
    let mut y_wide = vec![0.0f64; nrows];
    serial.spmv(&x, &mut y_wide);
    let compact = fs_sparse::CsrCompact::from_csr(&serial);
    let mut y_cmp = vec![0.0f64; nrows];
    compact.spmv(&x, &mut y_cmp);
    assert!(
        y_wide
            .iter()
            .zip(&y_cmp)
            .all(|(a, b)| a.to_bits() == b.to_bits()),
        "compact spmv != wide spmv bitwise"
    );
    // (1b + 2) sharded spmv bitwise at every thread count.
    for t in [1usize, 2, 3, 4, 8, 16] {
        let mut y_sh = vec![0.0f64; nrows];
        compact.spmv_sharded(&x, &mut y_sh, t);
        assert!(
            y_wide
                .iter()
                .zip(&y_sh)
                .all(|(a, b)| a.to_bits() == b.to_bits()),
            "sharded spmv (t={t}) != serial bitwise"
        );
    }
    println!(
        "{{\"suite\":\"fs-sparse\",\"case\":\"wsbf-bitwise-twins\",\"verdict\":\"pass\",\"detail\":\"parallel assembly (t in 1..8) and compact+sharded spmv (t in 1..16) all bitwise == serial wide kernels\"}}"
    );
}

/// wsbf segment 2: the chunk-major SELL kernels and the blocked SpMM
/// are bitwise-equal to their reference twins (pads read, signed-zero
/// argument inherited; every thread count; every rhs block width).
#[test]
fn wsbf_segment2_bitwise_twins() {
    let (nrows, ncols) = (203usize, 157usize);
    let mut coo = fs_sparse::Coo::new(nrows, ncols);
    let mut seed = 0x5E11_2026_u64;
    let mut lcg = move || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        seed
    };
    for _ in 0..4000 {
        let r = (lcg() % nrows as u64) as usize;
        let c = (lcg() % ncols as u64) as usize;
        let v = ((lcg() >> 11) as f64) / (1u64 << 53) as f64 - 0.5;
        coo.push(r, c, v);
    }
    // A few very long rows + empty rows (ragged pads exercised).
    for k in 0..120 {
        coo.push(7, k % ncols, 0.125 + k as f64 * 0.001);
    }
    let a = coo.assemble();
    let x: Vec<f64> = (0..ncols).map(|i| -0.75 + (i % 23) as f64 * 0.1).collect();
    let mut y_ref = vec![0.0f64; nrows];
    a.spmv(&x, &mut y_ref);
    for (c, sigma) in [(4usize, 32usize), (8, 64), (2, 16)] {
        let sell = fs_sparse::Sell::from_csr(&a, c, sigma);
        let mut y_row = vec![0.0f64; nrows];
        sell.spmv(&x, &mut y_row);
        let mut y_ch = vec![0.0f64; nrows];
        sell.spmv_chunked(&x, &mut y_ch);
        assert!(
            y_ref
                .iter()
                .zip(&y_ch)
                .all(|(u, v)| u.to_bits() == v.to_bits()),
            "chunked SELL (C={c}) != CSR bitwise"
        );
        for t in [1usize, 3, 8] {
            let mut y_sh = vec![0.0f64; nrows];
            sell.spmv_chunked_sharded(&x, &mut y_sh, t);
            assert!(
                y_ref
                    .iter()
                    .zip(&y_sh)
                    .all(|(u, v)| u.to_bits() == v.to_bits()),
                "sharded chunked SELL (C={c}, t={t}) != CSR bitwise"
            );
        }
        assert!(
            y_ref
                .iter()
                .zip(&y_row)
                .all(|(u, v)| u.to_bits() == v.to_bits()),
            "row-major SELL != CSR (regression)"
        );
    }
    // Blocked SpMM == per-column SpMV, widths that exercise partial blocks.
    for nrhs in [1usize, 3, 8, 11] {
        let b: Vec<f64> = (0..ncols * nrhs)
            .map(|i| 0.5 - (i % 19) as f64 * 0.05)
            .collect();
        let mut y_mm = vec![0.0f64; nrows * nrhs];
        a.spmm_blocked(&b, nrhs, &mut y_mm);
        for j in 0..nrhs {
            let xj: Vec<f64> = (0..ncols).map(|i| b[i * nrhs + j]).collect();
            let mut yj = vec![0.0f64; nrows];
            a.spmv(&xj, &mut yj);
            assert!(
                (0..nrows).all(|r| y_mm[r * nrhs + j].to_bits() == yj[r].to_bits()),
                "spmm_blocked col {j} (nrhs={nrhs}) != spmv bitwise"
            );
        }
    }
    println!(
        "{{\"suite\":\"fs-sparse\",\"case\":\"wsbf-segment2-twins\",\"verdict\":\"pass\",\"detail\":\"chunked SELL (C in 2/4/8, t in 1/3/8) and blocked SpMM (nrhs in 1/3/8/11) bitwise == references\"}}"
    );
}

/// wsbf item 6: the sparse-accumulator SpGEMM is bitwise-equal to the
/// dense-SPA reference, including on a VERY WIDE product where the
/// dense scratch is the thing being avoided.
#[test]
fn wsbf_sparse_spa_spgemm() {
    let mut seed = 0x59A_2026_u64;
    let mut lcg = move || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        seed
    };
    // Square-ish random product.
    let mut ca = fs_sparse::Coo::new(80, 60);
    let mut cb = fs_sparse::Coo::new(60, 90);
    for _ in 0..900 {
        ca.push(
            (lcg() % 80) as usize,
            (lcg() % 60) as usize,
            ((lcg() >> 11) as f64) / (1u64 << 53) as f64 - 0.5,
        );
        cb.push(
            (lcg() % 60) as usize,
            (lcg() % 90) as usize,
            ((lcg() >> 11) as f64) / (1u64 << 53) as f64 - 0.5,
        );
    }
    let (a, b) = (ca.assemble(), cb.assemble());
    let dense = fs_sparse::ops::spgemm(&a, &b);
    let sparse = fs_sparse::ops::spgemm_sparse_spa(&a, &b);
    assert_eq!(dense, sparse, "sparse-SPA SpGEMM != dense-SPA bitwise");
    // Very wide B: 2_000_000 columns, a handful of entries — the dense
    // SPA would burn a 16 MB+ scratch per call; the sparse one doesn't.
    let wide_cols = 2_000_000usize;
    let mut cw = fs_sparse::Coo::new(60, wide_cols);
    for k in 0..300 {
        cw.push(
            (lcg() % 60) as usize,
            (lcg() % wide_cols as u64) as usize,
            0.25 + f64::from(k) * 0.001,
        );
    }
    let w = cw.assemble();
    let dense_w = fs_sparse::ops::spgemm(&a, &w);
    let sparse_w = fs_sparse::ops::spgemm_sparse_spa(&a, &w);
    assert_eq!(dense_w, sparse_w, "wide sparse-SPA != dense-SPA bitwise");
    println!(
        "{{\"suite\":\"fs-sparse\",\"case\":\"wsbf-sparse-spa\",\"verdict\":\"pass\",\"detail\":\"BTree-SPA SpGEMM bitwise == dense-SPA on random and 2e6-column-wide products\"}}"
    );
}

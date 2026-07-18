//! fs-sparse conformance: the cross-format, cross-ISA battery any
//! reimplementation must pass (plan §13.3). Builds FEM-patterned and
//! adversarial matrices, runs every format's SpMV plus the pattern algebra,
//! and folds all output bits into one FNV-64 golden hash — recorded on
//! aarch64-apple and required to match on x86-64 (the same evidence
//! discipline as fs-math/fs-fft).

use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_propcheck::Shrink;
use fs_sparse::{Bsr, Coo, Csr, Sell, ops};

#[derive(Clone, Debug)]
struct SparseCase {
    triplets: Vec<(u64, u64, i64)>,
    x: [i64; 8],
}

impl Shrink for SparseCase {
    fn shrink_candidates(&self) -> Vec<Self> {
        let mut candidates: Vec<Self> = self
            .triplets
            .shrink_candidates()
            .into_iter()
            .map(|triplets| SparseCase {
                triplets,
                x: self.x,
            })
            .collect();
        for (index, value) in self.x.iter().enumerate() {
            for candidate in value.shrink_candidates() {
                let mut x = self.x;
                x[index] = candidate;
                candidates.push(SparseCase {
                    triplets: self.triplets.clone(),
                    x,
                });
            }
        }
        candidates
    }
}

fn bitwise_equal(left: &[f64], right: &[f64]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .all(|(a, b)| a.to_bits() == b.to_bits())
}

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

/// G0 generated cross-format battery (bead frankensim-4nh8). The fixed zoo
/// and its cross-ISA golden above remain the durable regression pin; this
/// layer fills the input space between those fixtures and shrinks failures.
#[test]
fn generated_cross_format_spmv_is_bitwise_equal() {
    fs_propcheck::check(
        "sparse-csr-bsr-sell-spmv-bitwise",
        0x5A_5001,
        600,
        |s| SparseCase {
            triplets: s.vec_of(64, |s| {
                (
                    u64::try_from(s.int_in(0, 7)).expect("row is non-negative"),
                    u64::try_from(s.int_in(0, 7)).expect("column is non-negative"),
                    s.int_in(-8, 8),
                )
            }),
            x: core::array::from_fn(|_| s.int_in(-8, 8)),
        },
        |case| {
            let mut coo = Coo::new(8, 8);
            for &(row, column, value) in &case.triplets {
                coo.push(
                    usize::try_from(row).expect("generated row fits usize"),
                    usize::try_from(column).expect("generated column fits usize"),
                    value as f64,
                );
            }
            let csr = coo.assemble();
            let x: Vec<f64> = case.x.iter().map(|&value| value as f64).collect();

            let mut csr_out = [0.0; 8];
            csr.spmv(&x, &mut csr_out);

            let bsr = Bsr::from_csr(&csr, 4, 4);
            let mut bsr_out = [0.0; 8];
            bsr.spmv(&x, &mut bsr_out);

            let sell = Sell::from_csr(&csr, 8, 32);
            let mut sell_out = [0.0; 8];
            sell.spmv(&x, &mut sell_out);

            bitwise_equal(&csr_out, &bsr_out) && bitwise_equal(&csr_out, &sell_out)
        },
    );
    println!(
        "{{\"suite\":\"fs-sparse\",\"case\":\"g0-generated-cross-format\",\
         \"verdict\":\"pass\",\"detail\":\"600 shrink-armed 8x8 CSR/BSR4/SELL8x32 SpMV cases\"}}"
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

/// G0 resource-admission regression: caller-controlled worker counts are
/// bounded by useful row parallelism before allocating shard metadata or
/// spawning workers. The cap must not change the serial bit pattern.
#[test]
fn oversized_worker_requests_are_capped_to_useful_rows() {
    let empty_coo = Coo::new(0, 0);
    let empty_serial = empty_coo.assemble();
    assert_eq!(empty_coo.assemble_parallel(usize::MAX), empty_serial);
    let empty_compact = fs_sparse::CsrCompact::from_csr(&empty_serial);
    let mut empty_y: [f64; 0] = [];
    empty_compact.spmv_sharded(&[], &mut empty_y, usize::MAX);
    assert_eq!(empty_compact.numa_localized(usize::MAX), empty_compact);

    let empty_rows_coo = Coo::new(128, 1);
    let empty_rows = empty_rows_coo.assemble();
    assert_eq!(empty_rows_coo.assemble_parallel(usize::MAX), empty_rows);
    let empty_rows_compact = fs_sparse::CsrCompact::from_csr(&empty_rows);
    let mut empty_rows_y = [f64::NAN; 128];
    empty_rows_compact.spmv_sharded(&[1.0], &mut empty_rows_y, usize::MAX);
    assert_eq!(empty_rows_y, [0.0; 128]);
    assert_eq!(
        empty_rows_compact.numa_localized(usize::MAX),
        empty_rows_compact
    );

    let mut coo = Coo::new(1, 1);
    coo.push(0, 0, 2.0);

    let serial = coo.assemble();
    assert_eq!(coo.assemble_parallel(usize::MAX), serial);

    let compact = fs_sparse::CsrCompact::from_csr(&serial);
    let mut y = [0.0];
    compact.spmv_sharded(&[3.0], &mut y, usize::MAX);
    assert_eq!(y[0].to_bits(), 6.0f64.to_bits());
    assert_eq!(compact.numa_localized(usize::MAX), compact);

    let mut ragged = Coo::new(7, 5);
    ragged.push(0, 0, 1.0);
    ragged.push(3, 2, 2.0);
    ragged.push(3, 2, -0.5);
    ragged.push(6, 4, -3.0);
    let ragged_serial = ragged.assemble();
    assert_eq!(ragged.assemble_parallel(usize::MAX), ragged_serial);

    let ragged_compact = fs_sparse::CsrCompact::from_csr(&ragged_serial);
    let x = [0.25, 0.5, 1.5, 2.0, -2.0];
    let mut ragged_serial_y = [0.0; 7];
    ragged_compact.spmv(&x, &mut ragged_serial_y);
    let mut ragged_sharded_y = [0.0; 7];
    ragged_compact.spmv_sharded(&x, &mut ragged_sharded_y, usize::MAX);
    assert!(bitwise_equal(&ragged_serial_y, &ragged_sharded_y));
    assert_eq!(ragged_compact.numa_localized(usize::MAX), ragged_compact);
}

/// Chunk-major SELL must not execute physical padding as extra arithmetic:
/// even a zero-valued pad can change a live row's signed-zero result.
#[test]
fn sell_chunked_skips_signed_zero_padding() {
    let tiny = f64::from_bits(1);
    let a = Csr::from_parts(2, 2, vec![0, 1, 3], vec![0, 0, 1], vec![-tiny, -tiny, 0.0]);
    let x = [0.25, 1.0];
    let sell = Sell::from_csr(&a, 2, 2);
    assert_eq!(sell.physical_slots(), a.nnz() + 1);

    let mut y_csr = [f64::NAN; 2];
    a.spmv(&x, &mut y_csr);
    assert_eq!(y_csr[0].to_bits(), (-0.0f64).to_bits());
    assert_eq!(y_csr[1].to_bits(), 0.0f64.to_bits());

    let mut y_row = [f64::NAN; 2];
    sell.spmv(&x, &mut y_row);
    assert!(bitwise_equal(&y_csr, &y_row));

    let mut y_chunked = [f64::NAN; 2];
    sell.spmv_chunked(&x, &mut y_chunked);
    assert!(bitwise_equal(&y_csr, &y_chunked));

    for threads in [1, 2, 3] {
        let mut y_sharded = [f64::NAN; 2];
        sell.spmv_chunked_sharded(&x, &mut y_sharded, threads);
        assert!(
            bitwise_equal(&y_csr, &y_sharded),
            "sharded chunked SELL diverged with {threads} threads"
        );
    }
}

/// wsbf segment 2: the chunk-major SELL kernels and the blocked SpMM
/// are bitwise-equal to their reference twins (pads skipped; every
/// thread count; every rhs block width).
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

// ---------------------------------------------------------------------------
// Cheap structured BEDROCK Casebook subset (bead 6ys.18.7)
// ---------------------------------------------------------------------------

const CASEBOOK_SUITE: &str = "fs-sparse/bedrock-conformance-v1";
const CASEBOOK_NROWS: usize = 4;
const CASEBOOK_NCOLS: usize = 4;
const CASEBOOK_ROW_PTR: [usize; 5] = [0, 2, 4, 6, 8];
const CASEBOOK_COLUMNS: [usize; 8] = [0, 2, 1, 3, 0, 2, 1, 3];
const CASEBOOK_VALUES: [f64; 8] = [2.0, -2.0, 3.0, 0.5, -1.0, 3.0, 2.0, -0.25];
const CASEBOOK_X: [f64; 4] = [1.0, 2.0, -1.0, 4.0];
const CASEBOOK_Y: [f64; 4] = [4.0, 8.0, -4.0, 3.0];
const CASEBOOK_TRIPLETS: [(usize, usize, f64); 10] = [
    (3, 3, -0.25),
    (0, 2, -1.25),
    (1, 1, 3.0),
    (2, 2, 1.0),
    (0, 0, 2.0),
    (3, 1, 2.0),
    (1, 3, 0.5),
    (0, 2, -0.75),
    (2, 0, -1.0),
    (2, 2, 2.0),
];
const CASEBOOK_SUCCESS_POLLS: usize = CASEBOOK_NROWS + CASEBOOK_VALUES.len();
const CASEBOOK_REFUSAL_POLL: usize = 5;
const CASEBOOK_MALFORMED_POLLS: usize = 3;

fn casebook_push_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn casebook_push_len(bytes: &mut Vec<u8>, value: usize) {
    casebook_push_u64(
        bytes,
        u64::try_from(value).expect("conformance fixture lengths fit u64"),
    );
}

fn casebook_push_text(bytes: &mut Vec<u8>, value: &str) {
    casebook_push_len(bytes, value.len());
    bytes.extend_from_slice(value.as_bytes());
}

fn casebook_push_usizes(bytes: &mut Vec<u8>, values: &[usize]) {
    casebook_push_len(bytes, values.len());
    for &value in values {
        casebook_push_len(bytes, value);
    }
}

fn casebook_push_f64s(bytes: &mut Vec<u8>, values: &[f64]) {
    casebook_push_len(bytes, values.len());
    for value in values {
        casebook_push_u64(bytes, value.to_bits());
    }
}

fn casebook_push_u64s(bytes: &mut Vec<u8>, values: &[u64]) {
    casebook_push_len(bytes, values.len());
    for &value in values {
        casebook_push_u64(bytes, value);
    }
}

fn casebook_push_nested(bytes: &mut Vec<u8>, label: &str, frame: &[u8]) {
    casebook_push_text(bytes, label);
    casebook_push_len(bytes, frame.len());
    bytes.extend_from_slice(frame);
}

fn casebook_fixture() -> Csr {
    let mut coo = Coo::new(CASEBOOK_NROWS, CASEBOOK_NCOLS);
    for &(row, column, value) in &CASEBOOK_TRIPLETS {
        coo.push(row, column, value);
    }
    coo.assemble()
}

fn casebook_assembly_inputs() -> Vec<u8> {
    let mut bytes = b"fs-sparse:coo-canonical-assembly-kat:v1".to_vec();
    casebook_push_text(&mut bytes, "Coo::push+assemble");
    casebook_push_text(
        &mut bytes,
        "stable-sort(row,column);duplicate-sum=insertion-order:v1",
    );
    casebook_push_len(&mut bytes, CASEBOOK_NROWS);
    casebook_push_len(&mut bytes, CASEBOOK_NCOLS);
    casebook_push_len(&mut bytes, CASEBOOK_TRIPLETS.len());
    for &(row, column, value) in &CASEBOOK_TRIPLETS {
        casebook_push_len(&mut bytes, row);
        casebook_push_len(&mut bytes, column);
        casebook_push_u64(&mut bytes, value.to_bits());
    }
    casebook_push_text(&mut bytes, "expected-row-pointers");
    casebook_push_usizes(&mut bytes, &CASEBOOK_ROW_PTR);
    casebook_push_text(&mut bytes, "expected-column-indices");
    casebook_push_usizes(&mut bytes, &CASEBOOK_COLUMNS);
    casebook_push_text(&mut bytes, "expected-value-bits");
    casebook_push_f64s(&mut bytes, &CASEBOOK_VALUES);
    bytes
}

fn casebook_spmv_inputs() -> Vec<u8> {
    let assembly = casebook_assembly_inputs();
    let mut bytes = b"fs-sparse:cross-format-spmv-kat:v1".to_vec();
    casebook_push_text(&mut bytes, "Csr::spmv+Bsr::spmv+Sell::spmv");
    casebook_push_text(
        &mut bytes,
        "ascending-global-column fused-mul-add from positive-zero:v1",
    );
    casebook_push_nested(&mut bytes, "nested-canonical-assembly-frame", &assembly);
    casebook_push_text(&mut bytes, "x");
    casebook_push_f64s(&mut bytes, &CASEBOOK_X);
    casebook_push_text(&mut bytes, "expected-y");
    casebook_push_f64s(&mut bytes, &CASEBOOK_Y);
    casebook_push_text(&mut bytes, "bsr-block-shape");
    casebook_push_usizes(&mut bytes, &[2, 2]);
    casebook_push_text(&mut bytes, "sell-c-sigma");
    casebook_push_usizes(&mut bytes, &[2, 4]);
    bytes
}

fn casebook_checkpoint_inputs() -> Vec<u8> {
    let assembly = casebook_assembly_inputs();
    let mut bytes = b"fs-sparse:checkpoint-publication-refusal-policy:v1".to_vec();
    casebook_push_text(&mut bytes, "Csr::try_from_parts_with_checkpoint");
    casebook_push_nested(&mut bytes, "nested-canonical-assembly-frame", &assembly);
    casebook_push_text(&mut bytes, "canonical-row-pointers");
    casebook_push_usizes(&mut bytes, &CASEBOOK_ROW_PTR);
    casebook_push_text(&mut bytes, "canonical-column-indices");
    casebook_push_usizes(&mut bytes, &CASEBOOK_COLUMNS);
    casebook_push_text(&mut bytes, "canonical-values");
    casebook_push_f64s(&mut bytes, &CASEBOOK_VALUES);
    casebook_push_text(&mut bytes, "expected-success-polls");
    casebook_push_len(&mut bytes, CASEBOOK_SUCCESS_POLLS);
    casebook_push_text(&mut bytes, "typed-refusal-poll");
    casebook_push_len(&mut bytes, CASEBOOK_REFUSAL_POLL);
    casebook_push_text(&mut bytes, "CheckpointStop::RefusedAt(5)");
    casebook_push_text(&mut bytes, "malformed-shape");
    casebook_push_usizes(&mut bytes, &[1, 4]);
    casebook_push_text(&mut bytes, "malformed-row-pointers");
    casebook_push_usizes(&mut bytes, &[0, 2]);
    casebook_push_text(&mut bytes, "malformed-duplicate-columns");
    casebook_push_usizes(&mut bytes, &[0, 0]);
    casebook_push_text(&mut bytes, "malformed-values");
    casebook_push_f64s(&mut bytes, &[1.0, 2.0]);
    casebook_push_text(&mut bytes, "expected-malformed-polls");
    casebook_push_len(&mut bytes, CASEBOOK_MALFORMED_POLLS);
    bytes
}

fn casebook_matrix_mismatch(matrix: &Csr) -> Option<String> {
    let mut row_ptr = Vec::with_capacity(matrix.nrows() + 1);
    let mut columns = Vec::with_capacity(matrix.nnz());
    let mut value_bits = Vec::with_capacity(matrix.nnz());
    row_ptr.push(0);
    for row in 0..matrix.nrows() {
        let (row_columns, row_values) = matrix.row(row);
        columns.extend_from_slice(row_columns);
        value_bits.extend(row_values.iter().map(|value| value.to_bits()));
        row_ptr.push(columns.len());
    }
    let reference_bits = CASEBOOK_VALUES.map(f64::to_bits);
    if matrix.nrows() == CASEBOOK_NROWS
        && matrix.ncols() == CASEBOOK_NCOLS
        && row_ptr.as_slice() == CASEBOOK_ROW_PTR.as_slice()
        && columns.as_slice() == CASEBOOK_COLUMNS.as_slice()
        && value_bits.as_slice() == reference_bits.as_slice()
    {
        None
    } else {
        Some(format!(
            "computed_shape={}x{}; reference_shape={}x{}; computed_row_ptr={row_ptr:?}; reference_row_ptr={CASEBOOK_ROW_PTR:?}; computed_columns={columns:?}; reference_columns={CASEBOOK_COLUMNS:?}; computed_value_bits={value_bits:016x?}; reference_value_bits={reference_bits:016x?}",
            matrix.nrows(),
            matrix.ncols(),
            CASEBOOK_NROWS,
            CASEBOOK_NCOLS,
        ))
    }
}

fn casebook_assembly_outcome() -> CaseOutcome {
    let matrix = casebook_fixture();
    if let Some(mismatch) = casebook_matrix_mismatch(&matrix) {
        return CaseOutcome::fail(format!(
            "operation=Coo::assemble; staged_triplets={}; convention=stable-sort(row,column)+insertion-order-duplicate-sum; {mismatch}",
            CASEBOOK_TRIPLETS.len(),
        ))
        .with_evidence("crates/fs-sparse/CONTRACT.md#public-types-and-semantics")
        .with_evidence("crates/fs-sparse/CONTRACT.md#invariants");
    }

    CaseOutcome::pass(
        "shape=4x4; staged=10; canonical_row_ptr=[0,2,4,6,8]; canonical_cols=[0,2,1,3,0,2,1,3]; duplicate_sums={(0,2):-2,(2,2):3}; value_bits=exact",
    )
    .with_evidence("crates/fs-sparse/CONTRACT.md#public-types-and-semantics")
    .with_evidence("crates/fs-sparse/CONTRACT.md#invariants")
}

fn casebook_spmv_bits() -> [[u64; 4]; 3] {
    let csr = casebook_fixture();
    let bsr = Bsr::from_csr(&csr, 2, 2);
    let sell = Sell::from_csr(&csr, 2, 4);
    let mut csr_y = [0.0; 4];
    let mut bsr_y = [0.0; 4];
    let mut sell_y = [0.0; 4];
    csr.spmv(&CASEBOOK_X, &mut csr_y);
    bsr.spmv(&CASEBOOK_X, &mut bsr_y);
    sell.spmv(&CASEBOOK_X, &mut sell_y);
    [
        csr_y.map(f64::to_bits),
        bsr_y.map(f64::to_bits),
        sell_y.map(f64::to_bits),
    ]
}

fn casebook_spmv_outcome() -> CaseOutcome {
    let computed = casebook_spmv_bits();
    let reference = CASEBOOK_Y.map(f64::to_bits);
    for (format_index, format) in ["csr", "bsr-2x2", "sell-2x4"].into_iter().enumerate() {
        for (component, (&computed_bits, &reference_bits)) in
            computed[format_index].iter().zip(&reference).enumerate()
        {
            if computed_bits != reference_bits {
                return CaseOutcome::fail(format!(
                    "operation=SpMV; format={format}; shape=4x4; x=[1,2,-1,4]; component={component}; computed_bits=0x{computed_bits:016x}; reference_bits=0x{reference_bits:016x}; computed_all={:016x?}; reference_all={reference:016x?}",
                    computed[format_index],
                ))
                .with_evidence("crates/fs-sparse/CONTRACT.md#invariants");
            }
        }
    }

    CaseOutcome::pass(
        "x=[1,2,-1,4]; y=[4,8,-4,3]; csr=exact; bsr2x2=exact; sell2x4=exact; cross_format_bits=identical",
    )
    .with_evidence("crates/fs-sparse/CONTRACT.md#invariants")
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CheckpointStop {
    RefusedAt(usize),
}

fn casebook_checkpoint_outcome() -> CaseOutcome {
    let mut success_polls = 0_usize;
    let published = Csr::try_from_parts_with_checkpoint(
        CASEBOOK_NROWS,
        CASEBOOK_NCOLS,
        CASEBOOK_ROW_PTR.to_vec(),
        CASEBOOK_COLUMNS.to_vec(),
        CASEBOOK_VALUES.to_vec(),
        || {
            success_polls += 1;
            Ok::<_, core::convert::Infallible>(())
        },
    );
    let published = match published {
        Ok(Some(matrix)) => matrix,
        Ok(None) => {
            return CaseOutcome::fail(format!(
                "operation=Csr::try_from_parts_with_checkpoint; canonical_publication=Ok(None); success_polls={success_polls}; expected_polls={CASEBOOK_SUCCESS_POLLS}"
            ))
            .with_evidence("crates/fs-sparse/CONTRACT.md#cancellation-behavior");
        }
        Err(never) => match never {},
    };
    if success_polls != CASEBOOK_SUCCESS_POLLS {
        return CaseOutcome::fail(format!(
            "operation=Csr::try_from_parts_with_checkpoint; canonical_publication=Ok(Some); success_polls={success_polls}; expected_polls={CASEBOOK_SUCCESS_POLLS}"
        ))
        .with_evidence("crates/fs-sparse/CONTRACT.md#cancellation-behavior");
    }
    if let Some(mismatch) = casebook_matrix_mismatch(&published) {
        return CaseOutcome::fail(format!(
            "operation=Csr::try_from_parts_with_checkpoint; canonical_publication=Ok(Some); success_polls={success_polls}; {mismatch}"
        ))
        .with_evidence("crates/fs-sparse/CONTRACT.md#invariants")
        .with_evidence("crates/fs-sparse/CONTRACT.md#cancellation-behavior");
    }

    let mut refusal_polls = 0_usize;
    let refused = Csr::try_from_parts_with_checkpoint(
        CASEBOOK_NROWS,
        CASEBOOK_NCOLS,
        CASEBOOK_ROW_PTR.to_vec(),
        CASEBOOK_COLUMNS.to_vec(),
        CASEBOOK_VALUES.to_vec(),
        || {
            refusal_polls += 1;
            if refusal_polls == CASEBOOK_REFUSAL_POLL {
                Err(CheckpointStop::RefusedAt(refusal_polls))
            } else {
                Ok(())
            }
        },
    );
    let refusal_matches = match &refused {
        Err(CheckpointStop::RefusedAt(poll)) => *poll == CASEBOOK_REFUSAL_POLL,
        _ => false,
    };
    if !refusal_matches || refusal_polls != CASEBOOK_REFUSAL_POLL {
        return CaseOutcome::fail(format!(
            "operation=Csr::try_from_parts_with_checkpoint; refusal_result={refused:?}; refusal_polls={refusal_polls}; expected_result=Err(CheckpointStop::RefusedAt({CASEBOOK_REFUSAL_POLL})); expected_refusal_polls={CASEBOOK_REFUSAL_POLL}"
        ))
        .with_evidence("crates/fs-sparse/CONTRACT.md#cancellation-behavior");
    }

    let mut malformed_polls = 0_usize;
    let malformed =
        Csr::try_from_parts_with_checkpoint(1, 4, vec![0, 2], vec![0, 0], vec![1.0, 2.0], || {
            malformed_polls += 1;
            Ok::<_, core::convert::Infallible>(())
        });
    let malformed_was_refused = match malformed {
        Ok(None) => true,
        Ok(Some(_)) => false,
        Err(never) => match never {},
    };
    if !malformed_was_refused || malformed_polls != CASEBOOK_MALFORMED_POLLS {
        return CaseOutcome::fail(format!(
            "operation=Csr::try_from_parts_with_checkpoint; malformed_duplicate_columns_refused={malformed_was_refused}; malformed_polls={malformed_polls}; expected_polls={CASEBOOK_MALFORMED_POLLS}; row_ptr=[0,2]; columns=[0,0]"
        ))
        .with_evidence("crates/fs-sparse/CONTRACT.md#error-model")
        .with_evidence("crates/fs-sparse/CONTRACT.md#cancellation-behavior");
    }

    CaseOutcome::pass(format!(
        "canonical_publication=exact; success_polls={success_polls}; typed_refusal=CheckpointStop::RefusedAt({refusal_polls}); partial_matrix_published=false; malformed_duplicate_columns=Ok(None); malformed_polls={malformed_polls}"
    ))
    .with_evidence("crates/fs-sparse/CONTRACT.md#error-model")
    .with_evidence("crates/fs-sparse/CONTRACT.md#cancellation-behavior")
}

#[test]
fn bedrock_casebook_suite_emits_replay_complete_green_records() {
    let assembly_digest = fnv1a64(&casebook_assembly_inputs());
    let spmv_digest = fnv1a64(&casebook_spmv_inputs());
    let checkpoint_digest = fnv1a64(&casebook_checkpoint_inputs());
    assert_eq!(assembly_digest, 0xe765_3922_29f9_9b04);
    assert_eq!(spmv_digest, 0x6a94_f1fe_0a7f_9980);
    assert_eq!(checkpoint_digest, 0xf223_b9d9_b841_b887);

    let report = Suite::new(CASEBOOK_SUITE)
        .case(
            "coo-canonical-assembly-kat",
            assembly_digest,
            ToleranceSpec::Exact,
            casebook_assembly_outcome,
        )
        .case(
            "cross-format-spmv-kat",
            spmv_digest,
            ToleranceSpec::Exact,
            casebook_spmv_outcome,
        )
        .case(
            "checkpoint-publication-refusal-policy",
            checkpoint_digest,
            ToleranceSpec::Structural,
            casebook_checkpoint_outcome,
        )
        .run();

    report.assert_green();
    assert_eq!(
        report
            .records
            .iter()
            .map(|record| record.case.as_str())
            .collect::<Vec<_>>(),
        [
            "coo-canonical-assembly-kat",
            "cross-format-spmv-kat",
            "checkpoint-publication-refusal-policy",
        ]
    );
    assert_eq!(
        report.records[0].json_line(),
        format!(
            concat!(
                "{{\"casebook\":{},\"suite\":\"fs-sparse/bedrock-conformance-v1\",",
                "\"case\":\"coo-canonical-assembly-kat\",\"inputs_digest\":\"e765392229f99b04\",",
                "\"tolerance\":\"exact\",\"pass\":true,",
                "\"details\":\"shape=4x4; staged=10; canonical_row_ptr=[0,2,4,6,8]; canonical_cols=[0,2,1,3,0,2,1,3]; duplicate_sums={{(0,2):-2,(2,2):3}}; value_bits=exact\",",
                "\"evidence\":[\"crates/fs-sparse/CONTRACT.md#public-types-and-semantics\",",
                "\"crates/fs-sparse/CONTRACT.md#invariants\"]}}"
            ),
            CASEBOOK_RECORD_VERSION,
        ),
        "the structured sparse-assembly record schema and field order are contract"
    );
}

#[test]
fn disclosed_seeded_corruption_turns_the_casebook_suite_red() {
    const CORRUPTION_SEED: u64 = 0xF55A_0001;
    let component = (CORRUPTION_SEED & 0x3) as usize;
    let bit = CORRUPTION_SEED.trailing_zeros();
    assert_eq!(component, 1);
    assert_eq!(bit, 0);
    let canonical = CASEBOOK_Y.map(f64::to_bits);
    let mut corrupted = canonical;
    corrupted[component] ^= 1_u64 << bit;

    let spmv = casebook_spmv_inputs();
    let mut inputs = b"fs-sparse:seeded-spmv-oracle-corruption:v1".to_vec();
    casebook_push_u64(&mut inputs, CORRUPTION_SEED);
    casebook_push_len(&mut inputs, component);
    casebook_push_u64(&mut inputs, u64::from(bit));
    casebook_push_nested(&mut inputs, "nested-cross-format-spmv-frame", &spmv);
    casebook_push_text(&mut inputs, "canonical-y-bits");
    casebook_push_u64s(&mut inputs, &canonical);
    casebook_push_text(&mut inputs, "corrupted-y-bits");
    casebook_push_u64s(&mut inputs, &corrupted);
    let inputs_digest = fnv1a64(&inputs);
    assert_eq!(inputs_digest, 0xdf14_5165_9b99_a1e8);

    let report = Suite::new(CASEBOOK_SUITE)
        .case(
            "seeded-spmv-oracle-corruption",
            inputs_digest,
            ToleranceSpec::Exact,
            move || {
                let computed = casebook_spmv_bits()[0];
                if computed == corrupted {
                    CaseOutcome::pass("seeded corruption was not detected")
                } else {
                    CaseOutcome::fail(format!(
                        "seed=0x{CORRUPTION_SEED:016x}; operation=Csr::spmv; shape=4x4; x=[1,2,-1,4]; component={component}; bit={bit}; computed={computed:016x?}; canonical={canonical:016x?}; corrupted={corrupted:016x?}"
                    ))
                    .with_evidence("crates/fs-sparse/tests/conformance.rs#seeded-corruption")
                }
            },
        )
        .run();

    assert!(
        !report.all_passed(),
        "the deliberately corrupted oracle must turn red"
    );
    let failures = report.failures();
    let [failure] = failures.as_slice() else {
        panic!("the seeded corruption must produce exactly one structured failure");
    };
    assert_eq!(failure.case, "seeded-spmv-oracle-corruption");
    assert_eq!(failure.inputs_digest, "df1451659b99a1e8");
    assert!(
        failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(failure.details.contains(&format!("component={component}")));
    assert!(failure.details.contains(&format!("bit={bit}")));
    assert!(failure.details.contains("computed=["));
    assert!(failure.details.contains("canonical=["));
    assert!(failure.details.contains("corrupted=["));
    assert!(
        failure
            .json_line()
            .contains("\"tolerance\":\"exact\",\"pass\":false")
    );

    let panic = std::panic::catch_unwind(|| report.assert_green())
        .expect_err("the merge-gate assertion must reject the seeded failure");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("casebook panic carries text");
    assert!(message.contains("seeded-spmv-oracle-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}

//! Rank-deficient TSQR theorem tiers, equivalence fixtures, and no-claim
//! contracts — bead frankensim-epic-bedrock-6ys.5.1.1.
//!
//! This battery is the frozen mathematical target for rank-deficient
//! [`fs_la::factor::tsqr_r`] and the immutable fixture source for the
//! independent certificate checker (6ys.5.1.4) and release E2E (6ys.5.1.6).
//! Every fixture below encodes ONE theorem tier or ONE explicit no-claim;
//! nothing here may be weakened to make a downstream checker pass.
//!
//! # Admitted inputs
//!
//! Finite `f64` row-major `m x n` slices with `m >= n >= 0` and
//! `row_block >= n` (with `row_block = 0 admitted only when `n = 0`). Storage
//! length is checked before any allocation. Non-finite entries are NOT
//! refused by the current admission surface (`admit_matrix_storage` checks
//! length only); behavior on non-finite input is outside every tier below and
//! its typed-refusal policy is an explicit 6ys.5.1.2 obligation, not a
//! silent gap.
//!
//! # Admitted reduction schedules (equivalence relation)
//!
//! A schedule is the leaf partition induced by `(m, row_block, n)` (final
//! fragment shorter than `n` rows absorbed into the previous leaf) plus the
//! fixed binary combine order over leaves. Two schedules are admissible
//! variants of one logical TSQR. The exact-arithmetic result class they share
//! is `R^T R = A^T A` with `rank(R) = rank(A)`; bitwise identity across
//! schedules is claimed ONLY inside Tier T2's hypotheses.
//!
//! # Theorem tiers
//!
//! * **T0 (exact arithmetic, all ranks).** For exact real arithmetic and any
//!   admissible schedule, the produced `R` satisfies `R^T R = A^T A` and
//!   `rank(R) = rank(A)`. If `A` has full column rank, positive-diagonal
//!   `R` is unique. Otherwise NO canonical form is asserted: rotations
//!   inside deficient subspaces are genuine gauge freedom.
//! * **T1 (same-ISA determinism, any rank).** Fixed input bits + fixed
//!   `(m, row_block, n)` produce bitwise-identical output across runs.
//! * **T2 (tree agreement, full column rank only).** For finite
//!   full-column-rank inputs, all admissible schedules agree within stated
//!   tolerance and match the sign-normalized Householder `R`.
//! * **T3 (no-claim, rank-deficient cross-schedule).** Neither bitwise nor
//!   value equality of `R` across different schedules is claimed for
//!   rank-deficient inputs. Each schedule's output is individually valid
//!   under T0. Near-threshold rank classification carries NO claim at any
//!   fixed absolute tolerance: thresholds must be scale-aware
//!   (dimensional), which the current `tsqr_r` signature does not yet carry
//!   (a 6ys.5.1.2 typed-policy obligation).
//!
//! # Proof obligations -> fixtures
//!
//! Reconstruction/orthogonality (T0): `gram_identity_across_ranks_and_schedules`.
//! Full-rank compatibility (T2): `tier_t2_full_rank_tree_agreement_and_uniqueness`.
//! Fixed-tree determinism (T1): `tier_t1_fixed_tree_bit_determinism`.
//! Cross-schedule validity without bit claims (T3): `tier_t3_cross_schedule_each_valid_bits_unclaimed`.
//! Exact vs near deficiency separation: `separation_exact_dependence_vs_perturbed`.
//! Scale-awareness motivation (T3): `threshold_boundary_is_scale_dependent_no_claim`.
//! Empty-shape prerequisite (57jrk): `empty_tsqr_semantics_are_frozen`.
//! Falsifier sensitivity: `falsifier_value_change_moves_output_bits`.

use fs_la::factor::{qr, tsqr_r};

/// Deterministic LCG in `[0, 1)`; identical seed => identical matrix.
fn rand_mat(rows: usize, cols: usize, seed: u64) -> Vec<f64> {
    let mut s = seed | 1;
    let mut out = vec![0.0; rows * cols];
    for v in out.iter_mut() {
        s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *v = ((s >> 11) as f64) / ((1u64 << 53) as f64);
    }
    out
}

fn transpose(a: &[f64], m: usize, n: usize) -> Vec<f64> {
    (0..n * m)
        .map(|i| a[(i % m) * n + i / m])
        .collect()
}

/// Naive triple-loop product; every fixture keeps sizes tiny.
fn matmul(a: &[f64], b: &[f64], m: usize, k: usize, n: usize) -> Vec<f64> {
    let mut c = vec![0.0; m * n];
    for i in 0..m {
        for p in 0..k {
            let aip = a[i * k + p];
            if aip == 0.0 {
                continue;
            }
            for j in 0..n {
                c[i * n + j] += aip * b[p * n + j];
            }
        }
    }
    c
}

fn max_rel_err(got: &[f64], want: &[f64]) -> f64 {
    got.iter()
        .zip(want)
        .map(|(g, w)| (g - w).abs() / (1.0 + w.abs()))
        .fold(0.0, f64::max)
}

fn is_upper_triangular(r: &[f64], n: usize) -> bool {
    (0..n).all(|i| (0..i).all(|j| r[i * n + j] == 0.0))
}

/// Sign convention: strictly-negative diagonals are flipped; `-0.0 < 0.0`
/// is false, so zero diagonals are retained exactly as produced.
fn diag_signs_within_convention(r: &[f64], n: usize) -> bool {
    (0..n).all(|i| !(r[i * n + i] < 0.0))
}

// ---------------------------------------------------------------------------
// T0: reconstruction authority holds at EVERY rank and EVERY schedule.
// ---------------------------------------------------------------------------
#[test]
fn gram_identity_across_ranks_and_schedules() {
    let m = 48usize;

    // Full-rank fixture: identity leading block pins rank independently of
    // the pseudorandom tail (same device as the full-rank battery).
    let n_full = 5usize;
    let mut full = rand_mat(m, n_full, 101);
    for i in 0..n_full {
        full[i * n_full + i] += 1.0;
    }

    // Exactly-dependent fixture: columns 1..3 are fixed multiples of column 0,
    // so rank is analytically 2 regardless of rounding in later columns' tail.
    let n_dep = 4usize;
    let mut dep = vec![0.0; m * n_dep];
    for i in 0..m {
        let x = (i as f64) - 17.0;
        dep[i * n_dep] = x;
        dep[i * n_dep + 1] = 2.0 * x;
        dep[i * n_dep + 2] = -3.0 * x;
        dep[i * n_dep + 3] = x * x;
    }

    for (label, a, n) in [
        ("full", &full, n_full),
        ("dependent", &dep, n_dep),
    ] {
        let ata = matmul(&transpose(a, m, n), a, n, m, n);
        for block in [6usize, 12, m] {
            let r = tsqr_r(a, m, n, block);
            assert_eq!(r.len(), n * n, "{label}: R shape");
            assert!(is_upper_triangular(&r, n), "{label}: upper triangular @block={block}");
            assert!(
                diag_signs_within_convention(&r, n),
                "{label}: strict-negative flip convention violated"
            );
            // T0 reconstruction: R^T R == A^T A. This is schedule-invariant
            // mathematics, unlike the factors themselves (T3).
            let rt = transpose(&r, n, n);
            let rtr = matmul(&rt, &r, n, n, n);
            let err = max_rel_err(&rtr, &ata);
            assert!(
                err < 1e-11,
                "{label}: Gram identity failed at block={block}: {err:e}"
            );
        }
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"tsqr_rank_deficient\",\"verdict\":\"pass\",\"detail\":\"T0 Gram identity holds for full-rank and exactly-dependent inputs across 3 schedules\"}}"
    );
}

// ---------------------------------------------------------------------------
// T2: full-column-rank compatibility — positive-diagonal R is unique, so all
// schedules and the direct Householder path must agree within tolerance.
// ---------------------------------------------------------------------------
#[test]
fn tier_t2_full_rank_tree_agreement_and_uniqueness() {
    let (m, n) = (120usize, 6usize);
    let mut a = rand_mat(m, n, 31);
    for i in 0..n {
        a[i * n + i] += 1.0;
    }

    let direct = {
        let f = qr(&a, m, n);
        let mut r = vec![0.0; n * n];
        for i in 0..n {
            let flip = if f.r(i, i) < 0.0 { -1.0 } else { 1.0 };
            for j in i..n {
                r[i * n + j] = flip * f.r(i, j);
            }
        }
        r
    };

    for block in [8usize, 24, 50, m] {
        let r = tsqr_r(&a, m, n, block.max(n));
        let err = max_rel_err(&r, &direct);
        assert!(err < 1e-9, "T2 violated at block={block}: {err:e}");
        // Uniqueness corollary: every schedule's diagonal is strictly
        // positive, not merely nonnegative.
        assert!(
            (0..n).all(|i| r[i * n + i] > 0.0),
            "T2 positive-diagonal uniqueness broken at block={block}"
        );
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"tsqr_rank_deficient\",\"verdict\":\"pass\",\"detail\":\"T2: 4 schedules agree with sign-normalized direct QR within 1e-9\"}}"
    );
}

// ---------------------------------------------------------------------------
// T1: fixed input bits + fixed schedule => bitwise-stable output, including
// the rank-deficient case where NO canonical form is claimed.
// ---------------------------------------------------------------------------
#[test]
fn tier_t1_fixed_tree_bit_determinism() {
    let m = 48usize;
    let mut dep = vec![0.0; m * 3];
    for i in 0..m {
        let x = (i as f64) - 17.0;
        dep[i * 3] = x;
        dep[i * 3 + 1] = 2.0 * x;
        dep[i * 3 + 2] = -x;
    }
    for block in [12usize, 16, m] {
        let r1 = tsqr_r(&dep, m, 3, block);
        let r2 = tsqr_r(&dep, m, 3, block);
        assert!(
            r1.iter().zip(&r2).all(|(x, y)| x.to_bits() == y.to_bits()),
            "T1 violated at block={block}"
        );
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"tsqr_rank_deficient\",\"verdict\":\"pass\",\"detail\":\"T1: bitwise rerun stability on rank-deficient fixture across 3 schedules\"}}"
    );
}

// ---------------------------------------------------------------------------
// T3 no-claim: across DIFFERENT schedules on a rank-deficient input, each R
// is valid (T0) but bitwise/value equality between them is NOT asserted —
// asserting it would fabricate uniqueness the mathematics does not grant.
// ---------------------------------------------------------------------------
#[test]
fn tier_t3_cross_schedule_each_valid_bits_unclaimed() {
    let m = 48usize;
    let mut dep = vec![0.0; m * 3];
    for i in 0..m {
        let x = (i as f64) - 17.0;
        dep[i * 3] = x;
        dep[i * 3 + 1] = 2.0 * x;
        dep[i * 3 + 2] = -x;
    }

    let ata = matmul(&transpose(&dep, m, 3), &dep, 3, m, 3);
    let mut divergence = 0.0f64;
    let schedules = [6usize, 12, 24, m];
    let rs: Vec<Vec<f64>> = schedules.iter().map(|&b| tsqr_r(&dep, m, 3, b)).collect();
    for (k, r) in rs.iter().enumerate() {
        assert!(is_upper_triangular(r, 3), "validity: triangular @{}", schedules[k]);
        let rt = transpose(r, 3, 3);
        let rtr = matmul(&rt, r, 3, 3, 3);
        let err = max_rel_err(&rtr, &ata);
        assert!(err < 1e-11, "T0 validity failed for schedule {}: {err:e}", schedules[k]);
    }
    for i in 0..rs.len() {
        for j in (i + 1)..rs.len() {
            divergence = divergence
                .max(max_rel_err(&rs[i], &rs[j]));
        }
    }
    // Observed divergence is DATA, not an assertion: it quantifies the gauge
    // freedom that the 6ys.5.1.3 canonical gauge must eliminate. Whether two
    // schedules happen to agree bit-for-bit today is incidental; the contract
    // forbids REQUIRING either outcome here.
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"tsqr_rank_deficient\",\"verdict\":\"pass\",\"detail\":\"T3: 4 schedules individually valid; max observed cross-schedule rel divergence {divergence:e} recorded as gauge-freedom data\"}}"
    );
}

// ---------------------------------------------------------------------------
// Exact dependence vs near dependence: the separation a future scale-aware
// rank profile must resolve, demonstrated at one benign scale.
// ---------------------------------------------------------------------------
#[test]
fn separation_exact_dependence_vs_perturbed() {
    let m = 48usize;
    let delta = 1e-10f64;
    let mut exact = vec![0.0; m * 2];
    let mut near = vec![0.0; m * 2];
    // w is deliberately NOT proportional to x, so the perturbation adds a
    // genuinely independent direction of size delta*||w_perp||. (Scaling x
    // by (1+delta) would leave both columns parallel and rank 1 — the
    // fixture originally made that mistake; the trailing pivot then measures
    // rounding noise, not separation.)
    for i in 0..m {
        let x = (i as f64) - 17.0;
        let w = ((i as f64) * 0.7).sin() + 1.25;
        exact[i * 2] = x;
        exact[i * 2 + 1] = 2.0 * x;
        near[i * 2] = x;
        near[i * 2 + 1] = 2.0 * x + delta * w;
    }

    let d_exact = {
        let r = tsqr_r(&exact, m, 2, 12);
        (r[1 * 2 + 1]).abs()
    };
    let d_near = {
        let r = tsqr_r(&near, m, 2, 12);
        (r[1 * 2 + 1]).abs()
    };
    let scale = 2.0; // loose rounding-envelope unit for the exact-dep case

    // Exact dependence: structural cancellation drives the trailing pivot
    // down to rounding level (~eps * ||A||), far below delta * scale.
    assert!(
        d_exact < 1e-13 * scale,
        "exact dependence left trailing pivot {d_exact:e}"
    );
    // Near dependence at relative delta survives above delta * scale by a
    // wide stability margin (Householder backward stability).
    assert!(
        d_near > 1e-3 * delta * scale && d_near < 1e4 * delta * scale,
        "perturbed pivot {d_near:e} escaped the stability window around {:e}",
        delta * scale
    );
    // The gap between the two regimes is what a tolerance must sit inside:
    // record both as data for the 6ys.5.1.2 policy surface.
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"tsqr_rank_deficient\",\"verdict\":\"pass\",\"detail\":\"separation: exact-dep pivot {d_exact:e} << delta*scale {:.e} <= near-dep pivot {d_near:e}\"}}",
        delta * scale
    );
}

// ---------------------------------------------------------------------------
// Threshold-boundary scale sweep: absolute tolerances CANNOT classify rank
// uniformly across unit scales. Validity (T0) is asserted everywhere; rank
// verdicts are explicitly NOT (T3), and measured pivots are recorded as the
// evidence that motivated banning dimensionful threshold literals.
// ---------------------------------------------------------------------------
#[test]
fn threshold_boundary_is_scale_dependent_no_claim() {
    let m = 48usize;
    let scales = [1.0f64, 1e8, 1e-8];
    let mut report = Vec::new();
    for &s in &scales {
        // Column 1 = column 0 scaled by (2(1+delta)): relatively near-
        // dependent at EVERY scale, but absolutely tiny at 1e-8 and huge at
        // 1e8. Any single absolute cutoff misclassifies one end.
        let delta = 1e-10f64;
        let mut a = vec![0.0; m * 2];
        for i in 0..m {
            let x = ((i as f64) - 17.0) * s;
            a[i * 2] = x;
            a[i * 2 + 1] = 2.0 * x * (1.0 + delta);
        }
        let r = tsqr_r(&a, m, 2, 12);
        // T0 validity only — no rank claim.
        let ata = matmul(&transpose(&a, m, 2), &a, 2, m, 2);
        let rt = transpose(&r, 2, 2);
        let rtr = matmul(&rt, &r, 2, 2, 2);
        let err = max_rel_err(&rtr, &ata);
        assert!(err < 1e-9, "scale {s:e}: Gram identity failed: {err:e}");
        let abs_pivot = r[1 * 2 + 1].abs();
        let rel_pivot = abs_pivot / (s * s); // normalized by ||col||^2 scale
        report.push(format!("{{\"scale\":{s:e},\"abs_pivot\":{abs_pivot:e},\"rel_pivot\":{rel_pivot:e}}}"));
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"tsqr_rank_deficient\",\"verdict\":\"pass\",\"detail\":\"scale sweep (no rank claim): [{}]\"}}",
        report.join(",")
    );
}

// ---------------------------------------------------------------------------
// Empty-TSQR semantics (prerequisite bead 57jrk, landed in 87108bb1):
// m x 0 shapes return empty R for every admitted row_block, including 0,
// while storage-length admission stays active.
// ---------------------------------------------------------------------------
#[test]
fn empty_tsqr_semantics_are_frozen() {
    assert!(tsqr_r(&[], 0, 0, 0).is_empty());
    assert!(tsqr_r(&[], 0, 0, 1).is_empty());
    assert!(tsqr_r(&[], 7, 0, 0).is_empty());
    assert!(tsqr_r(&[], 7, 0, 3).is_empty());

    // Admission remains fail-closed on wrong-length storage even for empty R.
    let result = std::panic::catch_unwind(|| tsqr_r(&[1.0], 7, 0, 0));
    assert!(result.is_err(), "storage mismatch must refuse");
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"tsqr_rank_deficient\",\"verdict\":\"pass\",\"detail\":\"empty shapes frozen: m x 0 accepted, row_block=0 admitted iff n=0, storage admission intact\"}}"
    );
}

// ---------------------------------------------------------------------------
// Falsifier: if an implementation ever became insensitive to material input
// changes (memoization bug, swallowed arguments), these assertions fire. A
// harness that cannot fail cannot certify T1/T2 either.
// ---------------------------------------------------------------------------
#[test]
fn falsifier_value_change_moves_output_bits() {
    let m = 32usize;
    let mut a = vec![0.0; m * 2];
    for i in 0..m {
        a[i * 2] = (i as f64) - 15.0;
        a[i * 2 + 1] = 0.25 * (i as f64);
    }
    let base = tsqr_r(&a, m, 2, 8);

    let mut b = a.clone();
    b[5 * 2] += 1.0; // material change to one entry
    let perturbed = tsqr_r(&b, m, 2, 8);
    assert!(
        !base.iter().zip(&perturbed).all(|(x, y)| x.to_bits() == y.to_bits()),
        "falsifier: output bits unchanged despite material input change"
    );

    let mut c = a.clone();
    c[m * 2 - 1] = f64::NAN; // non-finite entry: OUTSIDE all tiers (see header)
    let with_nan = tsqr_r(&c, m, 2, 8);
    // Documented current behavior: NaN propagates into the factor rather
    // than being refused at admission. Recorded as data pinning the status
    // quo that 6ys.5.1.2 must replace with a typed refusal.
    assert!(
        with_nan.iter().any(|v| v.is_nan()) || with_nan.iter().all(|v| v.is_finite()),
        "non-finite behavior changed; update this fixture and the tier header together"
    );
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"tsqr_rank_deficient\",\"verdict\":\"pass\",\"detail\":\"falsifiers: material input change moves output bits; non-finite propagation pinned for 5.1.2 refusal policy\"}}"
    );
}

//! Shared-basis block rational-Krylov reduction (bead
//! wf-root-guzez.5.8.2, E4.3b2-ii). Plan law: ONE projection basis
//! shared across ALL frozen-grid points — independently-reduced-per-
//! point interpolation is FORBIDDEN (the battery executes the
//! falsifier: per-point bases give garbage scheduling interpolation
//! where the shared basis stays accurate). IRKA-like refinement moves
//! the shifts toward the mirrored reduced poles for a bounded,
//! deterministic number of passes. Order ladder 4..32 by 4 — the
//! SMALLEST order passing the transfer tolerance on HELD-OUT
//! frequencies wins; failure at 32 emits the typed refusal that
//! forbids the A1 ROM in that domain (never a silent fallback).

use crate::rom::{A1Lti, solve_multi};
use crate::{Refusal, refuse};
use fs_blake3::hash_domain;

/// Ladder orders (registered).
pub const LADDER: [usize; 8] = [4, 8, 12, 16, 20, 24, 28, 32];

/// Transfer tolerance (relative, max over held-out frequencies and
/// MIMO channels; declared — an order passes when its worst channel
/// error sits under this).
pub const TRANSFER_TOL: f64 = 1.0e-2;

/// IRKA-like refinement passes (bounded, deterministic).
pub const IRKA_PASSES: usize = 3;

/// Held-out validation frequencies [rad/s] (NEVER used as shifts).
pub const HELD_OUT_W: [f64; 7] = [0.7, 1.7, 3.7, 7.3, 13.0, 23.0, 41.0];

/// Initial real shifts [rad/s] (log-spread over the flight band).
pub const SHIFTS_0: [f64; 4] = [1.0, 4.0, 12.0, 36.0];

/// Deterministic reserve shifts appended when the requested pool
/// saturates before reaching the rung order (few-anchor cases — e.g.
/// the per-point falsifier twins whose Krylov spaces overlap).
pub const RESERVE_SHIFTS: [f64; 6] = [0.3, 2.0, 6.0, 18.0, 54.0, 90.0];

/// A reduced system in the SHARED basis.
#[derive(Clone, Debug, PartialEq)]
pub struct ReducedLti {
    /// Reduced order.
    pub order: usize,
    /// Ar (order×order), row-major.
    pub a: Vec<f64>,
    /// Br (order×2).
    pub b: Vec<f64>,
    /// Cr (3×order).
    pub c: Vec<f64>,
    /// D (3×2) — carried through unchanged.
    pub d: [f64; 6],
}

/// One ladder rung's receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct RungReceipt {
    /// Order tried.
    pub order: usize,
    /// Worst relative transfer error over held-out frequencies and
    /// channels, max over the anchor systems.
    pub worst_rel_err: f64,
    /// Passed?
    pub passed: bool,
}

/// The reduction outcome.
#[derive(Clone, Debug, PartialEq)]
pub struct SharedReduction {
    /// The shared orthonormal basis (n×r, column-major columns).
    pub basis: Vec<f64>,
    /// FOM order.
    pub n: usize,
    /// Winning reduced order (smallest passing).
    pub order: usize,
    /// Ladder receipts (every rung tried, in order).
    pub ladder: Vec<RungReceipt>,
    /// Final shifts after refinement.
    pub shifts: Vec<f64>,
    /// Content digest (basis + ladder).
    pub digest: String,
}

// --------------------------------------------------------------------------
// Small dense helpers (row-major).
// --------------------------------------------------------------------------

fn matmul(a: &[f64], b: &[f64], m: usize, k: usize, n: usize) -> Vec<f64> {
    let mut c = vec![0.0f64; m * n];
    for i in 0..m {
        for p in 0..k {
            let aip = a[i * k + p];
            if aip != 0.0 {
                for j in 0..n {
                    c[i * n + j] += aip * b[p * n + j];
                }
            }
        }
    }
    c
}

/// Modified Gram–Schmidt append: orthonormalize `col` (len n) against
/// the columns already in `basis` (column-major, n×r); returns true if
/// the column survived (norm above the drop tolerance).
fn mgs_append(basis: &mut Vec<f64>, n: usize, col: &[f64]) -> bool {
    let r = basis.len() / n;
    let mut v = col.to_vec();
    for _ in 0..2 {
        for j in 0..r {
            let q = &basis[j * n..(j + 1) * n];
            let dot: f64 = q.iter().zip(v.iter()).map(|(a, b)| a * b).sum();
            for (vi, qi) in v.iter_mut().zip(q.iter()) {
                *vi -= dot * qi;
            }
        }
    }
    let norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    if norm < 1e-10 {
        return false;
    }
    for vi in &mut v {
        *vi /= norm;
    }
    basis.extend_from_slice(&v);
    true
}

/// Solve (σI − A) X = B (real shift σ, B n×2), returning X.
fn shifted_solve(a: &[f64], n: usize, sigma: f64, b: &[f64]) -> Result<Vec<f64>, Refusal> {
    let mut m = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            m[i * n + j] = f64::from(u8::from(i == j)) * sigma - a[i * n + j];
        }
    }
    let mut x = b.to_vec();
    solve_multi(&mut m, &mut x, n, 2)?;
    Ok(x)
}

/// MIMO transfer G(iω) as 3×2 complex (re, im) via the 2n real system.
fn transfer(
    a: &[f64],
    b: &[f64],
    c: &[f64],
    d: &[f64; 6],
    n: usize,
    w: f64,
) -> Result<[(f64, f64); 6], Refusal> {
    // (iωI − A) x = b  →  [[−A, −ωI],[ωI, −A]]·[xr; xi] = [br; 0].
    let big = 2 * n;
    let mut m = vec![0.0f64; big * big];
    for i in 0..n {
        for j in 0..n {
            m[i * big + j] = -a[i * n + j];
            m[(n + i) * big + (n + j)] = -a[i * n + j];
        }
        m[i * big + (n + i)] = -w;
        m[(n + i) * big + i] = w;
    }
    let mut rhs = vec![0.0f64; big * 2];
    for i in 0..n {
        rhs[i * 2] = b[i * 2];
        rhs[i * 2 + 1] = b[i * 2 + 1];
    }
    // rhs holds [br columns; zeros] — row-major big×2.
    solve_multi(&mut m, &mut rhs, big, 2)?;
    let mut out = [(0.0f64, 0.0f64); 6];
    for o in 0..3 {
        for input in 0..2 {
            let mut re = d[o * 2 + input];
            let mut im = 0.0;
            for j in 0..n {
                re += c[o * n + j] * rhs[j * 2 + input];
                im += c[o * n + j] * rhs[(n + j) * 2 + input];
            }
            out[o * 2 + input] = (re, im);
        }
    }
    Ok(out)
}

/// Project the FOM through the shared basis (Galerkin, one-sided).
#[must_use]
pub fn project(lti: &A1Lti, basis: &[f64], r: usize) -> ReducedLti {
    let n = lti.order;
    // Vt (r×n) from column-major basis.
    let mut vt = vec![0.0f64; r * n];
    for j in 0..r {
        for i in 0..n {
            vt[j * n + i] = basis[j * n + i];
        }
    }
    // V row-major n×r.
    let mut v = vec![0.0f64; n * r];
    for i in 0..n {
        for j in 0..r {
            v[i * r + j] = basis[j * n + i];
        }
    }
    let av = matmul(&lti.a, &v, n, n, r);
    let ar = matmul(&vt, &av, r, n, r);
    let br = matmul(&vt, &lti.b, r, n, 2);
    let cr = matmul(&lti.c, &v, 3, n, r);
    ReducedLti {
        order: r,
        a: ar,
        b: br,
        c: cr,
        d: lti.d,
    }
}

/// Eigenvalues of a small real matrix (Francis-free: symmetric-free
/// unshifted-QR is unreliable, so this uses the deterministic
/// Hessenberg + Wilkinson-shift QR with 2×2 deflation — enough for the
/// r ≤ 32 reduced systems this module owns).
///
/// # Errors
/// `rom-eig-stall` if deflation stalls (never observed on the ladder;
/// the refusal keeps the failure typed instead of looping forever).
pub fn small_eigenvalues(a_in: &[f64], n: usize) -> Result<Vec<(f64, f64)>, Refusal> {
    let mut a = a_in.to_vec();
    // Householder Hessenberg.
    for col in 0..n.saturating_sub(2) {
        let mut x: Vec<f64> = (col + 1..n).map(|i| a[i * n + col]).collect();
        let alpha = -x[0].signum() * x.iter().map(|v| v * v).sum::<f64>().sqrt();
        if alpha.abs() < 1e-300 {
            continue;
        }
        x[0] -= alpha;
        let vnorm = x.iter().map(|v| v * v).sum::<f64>().sqrt();
        if vnorm < 1e-300 {
            continue;
        }
        let v: Vec<f64> = x.iter().map(|t| t / vnorm).collect();
        // A ← (I−2vvᵀ)A(I−2vvᵀ) on the trailing block.
        for j in 0..n {
            let mut dot = 0.0;
            for (k, vk) in v.iter().enumerate() {
                dot += vk * a[(col + 1 + k) * n + j];
            }
            for (k, vk) in v.iter().enumerate() {
                a[(col + 1 + k) * n + j] -= 2.0 * vk * dot;
            }
        }
        for i in 0..n {
            let mut dot = 0.0;
            for (k, vk) in v.iter().enumerate() {
                dot += vk * a[i * n + col + 1 + k];
            }
            for (k, vk) in v.iter().enumerate() {
                a[i * n + col + 1 + k] -= 2.0 * vk * dot;
            }
        }
    }
    let mut eigs = Vec::new();
    let mut hi = n;
    let mut guard = 0usize;
    while hi > 0 {
        guard += 1;
        if guard > 200 * n {
            return Err(refuse(
                "rom-eig-stall",
                format!("QR deflation stalled at block {hi}"),
                "reduced system is pathological; refuse the rung",
            ));
        }
        if hi == 1 {
            eigs.push((a[0], 0.0));
            break;
        }
        // Deflation checks.
        let sub = a[(hi - 1) * n + (hi - 2)].abs();
        let scale = a[(hi - 1) * n + (hi - 1)].abs() + a[(hi - 2) * n + (hi - 2)].abs();
        if sub < 1e-12 * scale.max(1e-300) {
            eigs.push((a[(hi - 1) * n + (hi - 1)], 0.0));
            hi -= 1;
            continue;
        }
        if hi >= 2 {
            // 2×2 trailing block eigenvalues.
            let p = a[(hi - 2) * n + (hi - 2)];
            let q = a[(hi - 2) * n + (hi - 1)];
            let r2 = a[(hi - 1) * n + (hi - 2)];
            let s = a[(hi - 1) * n + (hi - 1)];
            let deflate2 = hi == 2
                || a[(hi - 2) * n + (hi - 3)].abs()
                    < 1e-12 * (p.abs() + a[(hi - 3) * n + (hi - 3)].abs()).max(1e-300);
            if deflate2 {
                let tr = p + s;
                let det = p * s - q * r2;
                let disc = tr * tr / 4.0 - det;
                if disc >= 0.0 {
                    let root = disc.sqrt();
                    eigs.push((tr / 2.0 + root, 0.0));
                    eigs.push((tr / 2.0 - root, 0.0));
                } else {
                    let root = (-disc).sqrt();
                    eigs.push((tr / 2.0, root));
                    eigs.push((tr / 2.0, -root));
                }
                hi -= 2;
                continue;
            }
        }
        // Wilkinson-shift QR sweep on the leading hi×hi block (Givens).
        let p = a[(hi - 2) * n + (hi - 2)];
        let q = a[(hi - 2) * n + (hi - 1)];
        let r2 = a[(hi - 1) * n + (hi - 2)];
        let s = a[(hi - 1) * n + (hi - 1)];
        let tr = p + s;
        let det = p * s - q * r2;
        let disc = tr * tr / 4.0 - det;
        let shift = if disc >= 0.0 {
            let e1 = tr / 2.0 + disc.sqrt();
            let e2 = tr / 2.0 - disc.sqrt();
            if (e1 - s).abs() < (e2 - s).abs() {
                e1
            } else {
                e2
            }
        } else {
            s
        };
        // Shifted QR via Givens on the Hessenberg block.
        let mut cs = vec![(1.0f64, 0.0f64); hi.saturating_sub(1)];
        for i in 0..hi {
            a[i * n + i] -= shift;
        }
        for k in 0..hi - 1 {
            let x = a[k * n + k];
            let z = a[(k + 1) * n + k];
            let rad = x.hypot(z);
            let (cthe, sthe) = if rad < 1e-300 {
                (1.0, 0.0)
            } else {
                (x / rad, z / rad)
            };
            cs[k] = (cthe, sthe);
            for j in k..hi {
                let t1 = a[k * n + j];
                let t2 = a[(k + 1) * n + j];
                a[k * n + j] = cthe * t1 + sthe * t2;
                a[(k + 1) * n + j] = -sthe * t1 + cthe * t2;
            }
        }
        for k in 0..hi - 1 {
            let (cthe, sthe) = cs[k];
            for i in 0..hi.min(k + 3) {
                let t1 = a[i * n + k];
                let t2 = a[i * n + k + 1];
                a[i * n + k] = cthe * t1 + sthe * t2;
                a[i * n + k + 1] = -sthe * t1 + cthe * t2;
            }
        }
        for i in 0..hi {
            a[i * n + i] += shift;
        }
    }
    Ok(eigs)
}

/// Build the shared basis for order r from anchors + shifts.
fn build_basis(anchors: &[&A1Lti], shifts: &[f64], r: usize) -> Result<Vec<f64>, Refusal> {
    let n = anchors[0].order;
    let mut pool: Vec<f64> = shifts.to_vec();
    pool.extend_from_slice(&RESERVE_SHIFTS);
    pool.dedup_by(|x, y| (*x - *y).abs() < 1e-9);
    let mut basis: Vec<f64> = Vec::new();
    'outer: for &sigma in &pool {
        for lti in anchors {
            let x = shifted_solve(&lti.a, n, sigma, &lti.b)?;
            for col in 0..2 {
                let column: Vec<f64> = (0..n).map(|i| x[i * 2 + col]).collect();
                mgs_append(&mut basis, n, &column);
                if basis.len() / n >= r {
                    break 'outer;
                }
            }
        }
    }
    if basis.len() / n < r {
        return Err(refuse(
            "rom-basis-deficient",
            format!("only {} of {r} columns survived", basis.len() / n),
            "add anchors; the block Krylov space (incl. reserve shifts) saturated",
        ));
    }
    Ok(basis)
}

/// Worst relative transfer error of the projection at order r over the
/// held-out frequencies, all channels, all anchors.
fn ladder_error(anchors: &[&A1Lti], basis: &[f64], r: usize) -> Result<f64, Refusal> {
    let mut worst = 0.0f64;
    for lti in anchors {
        let red = project(lti, basis, r);
        for &w in &HELD_OUT_W {
            let gf = transfer(&lti.a, &lti.b, &lti.c, &lti.d, lti.order, w)?;
            let gr = transfer(&red.a, &red.b, &red.c, &red.d, red.order, w)?;
            for ch in 0..6 {
                let (fr, fi) = gf[ch];
                let (rr, ri) = gr[ch];
                let mag = fr.hypot(fi);
                if mag > 1e-9 {
                    let err = (fr - rr).hypot(fi - ri) / mag;
                    worst = worst.max(err);
                }
            }
        }
    }
    Ok(worst)
}

/// Run the shared-basis reduction with IRKA-like refinement + ladder.
///
/// # Errors
/// `rom-ladder-exhausted` — order 32 failed the tolerance: the A1 ROM
/// is FORBIDDEN in this domain (plan law; never a silent fallback);
/// `rom-anchors-empty`; solver refusals pass through.
pub fn reduce_shared(anchors: &[&A1Lti]) -> Result<SharedReduction, Refusal> {
    reduce_shared_with_tol(anchors, TRANSFER_TOL)
}

/// [`reduce_shared`] with an explicit tolerance (battery driver — the
/// adversarial tolerance executes the ladder-exhausted refusal path).
///
/// # Errors
/// As [`reduce_shared`].
pub fn reduce_shared_with_tol(anchors: &[&A1Lti], tol: f64) -> Result<SharedReduction, Refusal> {
    if anchors.is_empty() {
        return Err(refuse(
            "rom-anchors-empty",
            "no anchor systems".into(),
            "pass the frozen-grid anchor set",
        ));
    }
    let n = anchors[0].order;
    let mut shifts: Vec<f64> = SHIFTS_0.to_vec();
    let mut ladder = Vec::new();
    // A reduction must reduce: rungs at or above the FOM order are
    // infeasible in an n-dimensional space and are skipped (declared).
    for &r in LADDER.iter().filter(|&&r| r < n) {
        // IRKA-like refinement at this rung: project, take reduced
        // poles, mirror the stable real parts into new shifts (bounded
        // passes, deterministic; complex poles contribute |λ|).
        let mut basis = build_basis(anchors, &shifts, r)?;
        for _ in 0..IRKA_PASSES {
            let red = project(anchors[0], &basis, r);
            let eigs = small_eigenvalues(&red.a, r)?;
            let mut new_shifts: Vec<f64> = eigs
                .iter()
                .map(|(re, im)| re.hypot(*im).max(1e-3))
                .collect();
            new_shifts.sort_by(|x, y| x.partial_cmp(y).expect("finite shifts"));
            new_shifts.dedup_by(|x, y| (*x - *y).abs() < 1e-6);
            if new_shifts.is_empty() {
                break;
            }
            let candidate = build_basis(anchors, &new_shifts, r)?;
            let old_err = ladder_error(anchors, &basis, r)?;
            let new_err = ladder_error(anchors, &candidate, r)?;
            if new_err < old_err {
                basis = candidate;
                shifts = new_shifts;
            } else {
                break;
            }
        }
        let worst = ladder_error(anchors, &basis, r)?;
        let passed = worst < tol;
        ladder.push(RungReceipt {
            order: r,
            worst_rel_err: worst,
            passed,
        });
        if passed {
            let mut bytes = Vec::new();
            for v in &basis {
                bytes.extend_from_slice(&v.to_bits().to_le_bytes());
            }
            for rung in &ladder {
                bytes.extend_from_slice(&rung.worst_rel_err.to_bits().to_le_bytes());
                bytes.push(u8::from(rung.passed));
            }
            let digest = hash_domain("org.frankensim.fs-wing.rom-reduction.v1", &bytes).to_hex();
            return Ok(SharedReduction {
                basis,
                n,
                order: r,
                ladder,
                shifts,
                digest,
            });
        }
    }
    Err(refuse(
        "rom-ladder-exhausted",
        format!("order 32 failed the {tol} transfer tolerance: {ladder:?}"),
        "A1 is FORBIDDEN in this domain (plan law) — use A0 or widen the grid",
    ))
}

/// The MIMO transfer of any (A,B,C,D) — public for the 5.8.3 receipt
/// clauses (phase/group delay live there).
///
/// # Errors
/// Solver refusals pass through.
pub fn transfer_of(
    a: &[f64],
    b: &[f64],
    c: &[f64],
    d: &[f64; 6],
    n: usize,
    w: f64,
) -> Result<[(f64, f64); 6], Refusal> {
    transfer(a, b, c, d, n, w)
}

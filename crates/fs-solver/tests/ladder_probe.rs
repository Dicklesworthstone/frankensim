//! Bead x08j acceptance: pMG-CG iteration counts must be FLAT across
//! BOTH ladders with a FIXED smoothing degree — the p-independence
//! (and h-robustness) the vertex-patch Schwarz smoother buys.
//!
//! Gate shape: the bead asks max/min <= 1.5 across r = 2..6 at fixed m
//! and across m = 2..4 at fixed r. At m = 2 there is exactly ONE
//! vertex patch covering the whole interior, so Schwarz is an EXACT
//! solve and CG converges in ~2 iterations — the literal min would
//! reward a trivial configuration. The m-ladder therefore gates the
//! ratio over the NONTRIVIAL meshes (m >= 3, extended to m = 5 and a
//! m = 6 spot-check — the window-sharing bug this caught at m >= 5 is
//! documented in pmg.rs) and separately requires the trivial m = 2
//! point to be at least as fast as m = 3.

use fs_solver::{CgState, MaskedTensorOp, PMultigrid};

fn rand_vec(n: usize, seed: u32) -> Vec<f64> {
    let mut s = u64::from(seed)
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (0..n)
        .map(|_| {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 11) as f64) / (1u64 << 53) as f64 - 0.5
        })
        .collect()
}

fn pmg_iters(m: usize, r: usize) -> usize {
    let op = MaskedTensorOp::new(m, r);
    let mut b = rand_vec(
        op.space().ndof(),
        50 + u32::try_from(10 * m + r).expect("small"),
    );
    for (bi, &mk) in b.iter_mut().zip(op.mask()) {
        if !mk {
            *bi = 0.0;
        }
    }
    // FIXED degree 3 across the whole ladder — the entire point: no
    // r-dependent crutch.
    let pmg = PMultigrid::new(m, r, 3);
    let mut st = CgState::new(&op, &pmg, &b);
    let rep = st.run(&op, &pmg, 1e-10, 100);
    assert!(rep.converged, "pMG-CG failed at m={m} r={r}: {rep:?}");
    rep.iters
}

#[test]
fn x08j_flat_ladders() {
    // ORDER ladder at fixed m = 3.
    let r_counts: Vec<usize> = (2..=6).map(|r| pmg_iters(3, r)).collect();
    let (rmin, rmax) = (
        *r_counts.iter().min().expect("nonempty"),
        *r_counts.iter().max().expect("nonempty"),
    );
    println!(
        "{{\"suite\":\"fs-solver\",\"case\":\"x08j-r-ladder\",\"verdict\":\"info\",\"detail\":\"m=3 r=2..6 iters={r_counts:?}\"}}"
    );
    assert!(
        2 * rmax <= 3 * rmin,
        "order ladder not flat (max/min > 1.5): {r_counts:?}"
    );
    // MESH ladder at fixed r = 3 (m = 2 is the trivially-exact
    // single-patch case — gated separately; ratio over m >= 3).
    let m_counts: Vec<usize> = (2..=5).map(|m| pmg_iters(m, 3)).collect();
    let nontrivial = &m_counts[1..];
    let (mmin, mmax) = (
        *nontrivial.iter().min().expect("nonempty"),
        *nontrivial.iter().max().expect("nonempty"),
    );
    println!(
        "{{\"suite\":\"fs-solver\",\"case\":\"x08j-m-ladder\",\"verdict\":\"info\",\"detail\":\"r=3 m=2..5 iters={m_counts:?}\"}}"
    );
    assert!(
        2 * mmax <= 3 * mmin,
        "mesh ladder not flat over m>=3 (max/min > 1.5): {m_counts:?}"
    );
    assert!(
        m_counts[0] <= m_counts[1],
        "the single-patch m=2 case must not be SLOWER than m=3: {m_counts:?}"
    );
    // Spot-check the regime that exposed the window-sharing bug
    // (three same-signature vertices per axis).
    let spot = pmg_iters(6, 2);
    assert!(
        2 * spot.max(mmax) <= 3 * spot.min(mmin).min(rmin),
        "m=6 spot-check out of family: {spot} vs ladders {r_counts:?} {m_counts:?}"
    );
    println!(
        "{{\"suite\":\"fs-solver\",\"case\":\"x08j-flat-ladders\",\"verdict\":\"pass\",\"detail\":\"r-ladder {r_counts:?} (ratio<=1.5), m-ladder {m_counts:?} (ratio<=1.5 over m>=3, m=2 exact-patch faster), m=6 spot {spot}\"}}"
    );
}

//! Bead avuw battery: Stokes saddle-point solves over the FEEC tensor
//! pair Q_r³/Q_{r−1}^disc with PMINRES + blockdiag(p-MG, pressure
//! mass). Gates: agreement with a dense pinned-pressure LU reference,
//! divergence-freeness of the velocity to solver tolerance, flat-ish
//! iteration counts across BOTH mesh and order ladders, bitwise
//! resume, and a frozen golden hash (separate constant — the tfz.10
//! solver golden is untouched).

use fs_la::factor::lu;
use fs_solver::{LinearOp, PminresState, StokesBlockDiag, StokesOp, StokesSystem, norm2};

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-solver\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

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

/// Masked random momentum forcing, zero mass source, projected rhs.
fn stokes_rhs(sys: &StokesSystem, seed: u32) -> Vec<f64> {
    let nv = sys.nv;
    let mut rhs = vec![0.0f64; sys.n()];
    let f = rand_vec(3 * nv, seed);
    for comp in 0..3 {
        for i in 0..nv {
            if sys.vmask()[i] {
                rhs[comp * nv + i] = f[comp * nv + i];
            }
        }
    }
    rhs
}

fn solve(sys: &StokesSystem, seed: u32, tol: f64) -> (Vec<f64>, usize, bool, f64) {
    let op = StokesOp::new(sys);
    let pre = StokesBlockDiag::new(sys, 3);
    let rhs = stokes_rhs(sys, seed);
    let mut st = PminresState::new(&op, &pre, &rhs);
    let rep = st.run(&op, &pre, tol, 400);
    let mut ax = vec![0.0f64; sys.n()];
    op.apply(&st.x, &mut ax);
    let diff: Vec<f64> = rhs.iter().zip(&ax).map(|(b, ax)| b - ax).collect();
    let true_rel = norm2(&diff) / norm2(&rhs).max(f64::MIN_POSITIVE);
    (st.x, rep.iters, rep.converged, true_rel)
}

#[test]
fn stokes_matches_dense_lu() {
    // m=2, r=2: reduced system = interior velocity dofs + pressure
    // with ONE dof pinned (the classic nonsingular reference).
    let sys = StokesSystem::new(2, 2);
    let op = StokesOp::new(&sys);
    let nv = sys.nv;
    let n = sys.n();
    let rhs = stokes_rhs(&sys, 77);
    // Reduced index set: interior velocities, all pressures except one
    // cell-constant coefficient.
    let mut red = Vec::new();
    for comp in 0..3 {
        for i in 0..nv {
            if sys.vmask()[i] {
                red.push(comp * nv + i);
            }
        }
    }
    // Pin pressure dof 0 — a P0 mode, whose constraint row is a linear
    // combination of the other cells' P0 rows (pinning a HIGH mode was
    // MEASURED to drop an independent constraint: pressure off by 100%).
    for p in 1..sys.np {
        red.push(3 * nv + p);
    }
    let nr = red.len();
    // Dense reduced matrix via unit-vector applies of the UNPROJECTED
    // saddle blocks (the projection is a solver device, not physics):
    // rebuild columns from StokesOp minus the projection by applying to
    // unit vectors and reading reduced rows. The projection only
    // touches the pressure block; on the pinned reduced system the
    // difference is rank-one in the null direction, which pinning
    // eliminates — but to stay exact we assemble from B directly.
    let mut kd = vec![0.0f64; nr * nr];
    let mut x = vec![0.0f64; n];
    let mut y = vec![0.0f64; n];
    for (cj, &dj) in red.iter().enumerate() {
        x.fill(0.0);
        x[dj] = 1.0;
        // Unprojected apply: velocity rows as in StokesOp; pressure
        // rows = B u (no projection).
        op.apply(&x, &mut y);
        // Undo the projection on the pressure block by recomputing it
        // raw: y_p = B u where u = x's velocity part. StokesOp
        // projected it; recover raw via divergence_inf-style apply.
        // (For unit vectors this is cheap and exact.)
        let raw_p = {
            let mut bu = vec![0.0f64; sys.np];
            // B is not public; use op on a velocity-only vector and
            // add back the projected component via the null vector.
            // Simpler: the projection subtracts (e·Bu)e, so
            // y_p + (e·Bu)e = Bu. We recover e·Bu from the pinned dof?
            // Not available — instead assemble from divergence of the
            // unit vector: e·Bu is the total "mass flux", ZERO for
            // interior velocity unit vectors on the enclosed domain
            // (columns of B sum against the constant to the boundary
            // integral, which vanishes). So y_p IS raw for velocity
            // columns; for pressure columns B-part is zero anyway.
            bu.copy_from_slice(&y[3 * nv..]);
            bu
        };
        for (ri, &di) in red.iter().enumerate() {
            kd[ri * nr + cj] = if di < 3 * nv {
                y[di]
            } else {
                raw_p[di - 3 * nv]
            };
        }
    }
    let fact = lu(&kd, nr);
    assert!(fact.is_ok(), "pinned Stokes reference is nonsingular");
    let Ok(fact) = fact else {
        return;
    };
    let mut xr: Vec<f64> = red.iter().map(|&d| rhs[d]).collect();
    fact.solve(&mut xr);
    // Iterative solution.
    let (xs, iters, conv, true_rel) = solve(&sys, 77, 1e-11);
    assert!(conv, "PMINRES must converge on the m=2 r=2 fixture");
    assert!(
        true_rel < 1e-9,
        "PMINRES estimate must agree with the true residual: {true_rel:.3e}"
    );
    // Compare velocities (unique) and mean-zero pressures.
    let mut dev = 0.0f64;
    let mut scale = 0.0f64;
    for (ri, &d) in red.iter().enumerate() {
        if d < 3 * nv {
            dev = dev.max((xr[ri] - xs[d]).abs());
            scale = scale.max(xr[ri].abs());
        }
    }
    assert!(
        dev < 1e-7 * scale.max(1.0),
        "velocity deviates from dense LU: {dev:.3e} (scale {scale:.3e})"
    );
    // Pressures agree after mean-zeroing both.
    let mut plu = vec![0.0f64; sys.np];
    for (ri, &d) in red.iter().enumerate() {
        if d >= 3 * nv {
            plu[d - 3 * nv] = xr[ri];
        }
    }
    let mut pit = xs[3 * nv..].to_vec();
    sys.project_pressure(&mut plu);
    sys.project_pressure(&mut pit);
    let pdev = plu
        .iter()
        .zip(&pit)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    let pscale = plu.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
    assert!(
        pdev < 1e-6 * pscale.max(1.0),
        "pressure deviates from dense LU: {pdev:.3e} (scale {pscale:.3e})"
    );
    log(
        "stokes-lu-agreement",
        "pass",
        &format!(
            "m=2 r=2: velocity dev {dev:.2e}, mean-zero pressure dev {pdev:.2e}, {iters} PMINRES iters"
        ),
    );
}

#[test]
fn stokes_velocity_is_divergence_free() {
    let sys = StokesSystem::new(3, 2);
    let (x, iters, conv, true_rel) = solve(&sys, 91, 1e-11);
    assert!(conv, "PMINRES must converge");
    assert!(
        true_rel < 1e-8,
        "PMINRES estimate must agree with the true residual: {true_rel:.3e}"
    );
    let div = sys.divergence_inf(&x);
    let uscale = x[..3 * sys.nv]
        .iter()
        .map(|v| v.abs())
        .fold(0.0f64, f64::max);
    assert!(
        div <= 1e-8 * uscale.max(1.0),
        "velocity not divergence-free: |Bu| = {div:.3e} vs scale {uscale:.3e}"
    );
    log(
        "stokes-div-free",
        "pass",
        &format!("m=3 r=2: |Bu|_inf = {div:.2e} at velocity scale {uscale:.2e} ({iters} iters)"),
    );
}

#[test]
fn stokes_iteration_envelope() {
    // Both ladders at fixed tolerance and smoothing degree: counts
    // must sit in a flat envelope (the whole point of the block
    // preconditioner — h- and p-robustness inherited from p-MG and the
    // diagonal pressure mass).
    let mut table = Vec::new();
    for &(m, r) in &[(2usize, 2usize), (3, 2), (4, 2), (2, 3), (3, 3)] {
        let sys = StokesSystem::new(m, r);
        let seed = 50 + u32::try_from(10 * m + r).unwrap_or(0);
        let (_, iters, conv, true_rel) = solve(&sys, seed, 1e-10);
        assert!(conv, "PMINRES failed at m={m} r={r}");
        assert!(
            true_rel < 1e-7,
            "true residual too large at m={m} r={r}: {true_rel:.3e}"
        );
        table.push((m, r, iters));
    }
    let counts: Vec<usize> = table.iter().map(|t| t.2).collect();
    let max = counts.iter().copied().max().unwrap_or(0);
    assert!(max <= 120, "Stokes iteration envelope exceeded: {table:?}");
    log(
        "stokes-envelope",
        "pass",
        &format!("(m, r, iters) = {table:?}; max {max} <= 120"),
    );
}

#[test]
fn stokes_resume_is_bitwise() {
    let sys = StokesSystem::new(2, 2);
    let op = StokesOp::new(&sys);
    let pre = StokesBlockDiag::new(&sys, 3);
    let rhs = stokes_rhs(&sys, 13);
    let mut straight = PminresState::new(&op, &pre, &rhs);
    straight.run(&op, &pre, 1e-11, 200);
    for cut in [1usize, 7, 23] {
        let mut a = PminresState::new(&op, &pre, &rhs);
        a.run(&op, &pre, 1e-11, cut);
        a.run(&op, &pre, 1e-11, 200 - cut);
        assert_eq!(a.iters, straight.iters, "iters differ at cut {cut}");
        for (u, v) in a.x.iter().zip(&straight.x) {
            assert_eq!(u.to_bits(), v.to_bits(), "resume not bitwise at cut {cut}");
        }
    }
    log("stokes-resume", "pass", "cuts 1/7/23 bitwise == straight");
}

// Frozen at bead avuw (independent of the tfz.10 solver golden).
const STOKES_GOLDEN: u64 = 0x5754_3908_cb41_7281;

#[test]
fn stokes_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let sys = StokesSystem::new(2, 3);
    let (x, iters, conv, true_rel) = solve(&sys, 21, 1e-10);
    assert!(conv);
    assert!(true_rel < 1e-7, "true residual too large: {true_rel:.3e}");
    for v in x.iter().step_by(11) {
        feed(*v);
    }
    #[allow(clippy::cast_precision_loss)]
    feed(iters as f64);
    log("stokes-golden", "info", &format!("0x{acc:016x}"));
    assert_eq!(
        acc, STOKES_GOLDEN,
        "Stokes bits changed: 0x{acc:016x} vs 0x{STOKES_GOLDEN:016x} — bump only with semantic justification"
    );
}

//! fs-solver battery (tfz.10 slices 1–2): CG/MINRES/GMRES correctness
//! against direct solves, bitwise pause/resume (the P7 obligation, at
//! each solver's stated granularity), transposed solves through the
//! same machinery (adjoint readiness), structured stall diagnoses,
//! deterministic repeat runs (G5), and p-multigrid: matrix-free
//! V-cycles with exact-injection transfers whose CG iteration counts
//! stay flat across BOTH the order ladder and the mesh ladder while
//! identity-preconditioned counts grow. Golden hash at the end.

use fs_feec::{element_geometry, kuhn_cube};
use fs_rand::StreamKey;
use fs_solver::{
    CgState, CsrF32, CsrOp, GmresState, LinearOp, MaskedTensorOp, MinresState, PMultigrid,
    StallDiagnosis, dot, mixed_cg_refine, norm2,
};
use fs_sparse::precond::IdentityPrecond;

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-solver\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

fn rand_vec(n: usize, tile: u32) -> Vec<f64> {
    let mut s = StreamKey {
        seed: 21,
        kernel: 0x501E,
        tile,
    }
    .stream();
    (0..n).map(|_| 2.0f64.mul_add(s.next_f64(), -1.0)).collect()
}

/// The interior-reduced Poisson CSR from fs-feec (kuhn(m), P1).
fn poisson_csr(m: usize) -> (fs_sparse::Csr, usize) {
    let (complex, positions) = kuhn_cube(m);
    let geo = element_geometry(&complex, &positions);
    let k0 = fs_feec::stiffness(
        &fs_feec::incidence_to_csr(&complex.d0()),
        &fs_feec::mass_matrix(&complex, &geo, 1),
    );
    let interior: Vec<usize> = (0..positions.len())
        .filter(|&v| !fs_feec::on_unit_cube_boundary(positions[v]))
        .collect();
    let mut slot = vec![usize::MAX; positions.len()];
    for (i, &v) in interior.iter().enumerate() {
        slot[v] = i;
    }
    let mut red = fs_sparse::Coo::new(interior.len(), interior.len());
    for (i, &v) in interior.iter().enumerate() {
        let (cols, vals) = k0.row(v);
        for (&c, &val) in cols.iter().zip(vals) {
            if slot[c] != usize::MAX {
                red.push(i, slot[c], val);
            }
        }
    }
    let n = interior.len();
    (red.assemble(), n)
}

#[test]
fn cg_solves_poisson_and_matches_direct() {
    let (a, n) = poisson_csr(2);
    let dense = a.to_dense();
    let op = CsrOp::symmetric(a);
    let b = rand_vec(n, 1);
    let mut st = CgState::new(&op, &IdentityPrecond, &b);
    let rep = st.run(&op, &IdentityPrecond, 1e-12, 10_000);
    assert!(rep.converged, "CG failed: {rep:?}");
    // Direct reference.
    let lu = fs_la::factor::lu(&dense, n).expect("nonsingular");
    let mut x_ref = b.clone();
    lu.solve(&mut x_ref);
    let dev =
        st.x.iter()
            .zip(&x_ref)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
    let scale = x_ref.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
    assert!(dev < 1e-9 * scale.max(1.0), "CG vs LU deviation {dev:.3e}");
    assert!(!rep.history.is_empty(), "history must be populated");
    log(
        "cg-poisson",
        "pass",
        &format!("iters={} dev={dev:.2e}", rep.iters),
    );
}

#[test]
fn cg_resume_is_bitwise() {
    let (a, n) = poisson_csr(4);
    let op = CsrOp::symmetric(a);
    let b = rand_vec(n, 2);
    let mut straight = CgState::new(&op, &IdentityPrecond, &b);
    straight.run(&op, &IdentityPrecond, 1e-13, 400);
    for cut in [1usize, 7, 40] {
        let mut first = CgState::new(&op, &IdentityPrecond, &b);
        first.run(&op, &IdentityPrecond, 1e-13, cut);
        let mut resumed = first.clone(); // checkpoint = clone
        resumed.run(&op, &IdentityPrecond, 1e-13, 400 - cut);
        assert_eq!(
            resumed.iters, straight.iters,
            "iter count differs at cut {cut}"
        );
        for (x1, x2) in resumed.x.iter().zip(&straight.x) {
            assert_eq!(x1.to_bits(), x2.to_bits(), "x bits differ at cut {cut}");
        }
    }
    // G5: repeat run bitwise.
    let mut again = CgState::new(&op, &IdentityPrecond, &b);
    again.run(&op, &IdentityPrecond, 1e-13, 400);
    assert!(
        again
            .x
            .iter()
            .zip(&straight.x)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
    log("cg-resume", "pass", "3 cut points + repeat bitwise");
}

#[test]
fn minres_handles_symmetric_indefinite() {
    // Indefinite diagonal-perturbed Poisson: A − σI with σ inside the
    // spectrum (CG's rz would go negative; MINRES is the right tool).
    let (a, n) = poisson_csr(4);
    let mut coo = fs_sparse::Coo::new(n, n);
    for r in 0..n {
        let (cols, vals) = a.row(r);
        for (&c, &v) in cols.iter().zip(vals) {
            coo.push(r, c, v);
        }
        coo.push(r, r, -12.0); // shift into indefiniteness
    }
    let shifted = coo.assemble();
    let dense = shifted.to_dense();
    let op = CsrOp::symmetric(shifted);
    let b = rand_vec(n, 3);
    let mut st = MinresState::new(&op, &b);
    let rep = st.run(&op, 1e-11, 10_000);
    assert!(rep.converged, "MINRES failed: {rep:?}");
    // True residual cross-check (the |η| estimate must be honest).
    let mut ax = vec![0.0f64; n];
    op.apply(&st.x, &mut ax);
    let true_rel = norm2(
        &b.iter()
            .zip(&ax)
            .map(|(bi, ai)| bi - ai)
            .collect::<Vec<_>>(),
    ) / norm2(&b);
    assert!(
        true_rel < 1e-9,
        "MINRES estimate dishonest: true rel {true_rel:.3e}"
    );
    // Direct reference.
    let lu = fs_la::factor::lu(&dense, n).expect("nonsingular");
    let mut x_ref = b.clone();
    lu.solve(&mut x_ref);
    let dev =
        st.x.iter()
            .zip(&x_ref)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
    let scale = x_ref.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
    assert!(
        dev < 1e-7 * scale.max(1.0),
        "MINRES vs LU deviation {dev:.3e}"
    );
    log(
        "minres-indefinite",
        "pass",
        &format!("iters={} true_rel={true_rel:.2e}", rep.iters),
    );
}

#[test]
fn minres_resume_is_bitwise() {
    let (a, n) = poisson_csr(4);
    let op = CsrOp::symmetric(a);
    let b = rand_vec(n, 4);
    let mut straight = MinresState::new(&op, &b);
    straight.run(&op, 1e-12, 300);
    for cut in [1usize, 11, 60] {
        let mut first = MinresState::new(&op, &b);
        first.run(&op, 1e-12, cut);
        let mut resumed = first.clone();
        resumed.run(&op, 1e-12, 300 - cut);
        assert_eq!(resumed.iters, straight.iters, "iters differ at cut {cut}");
        for (x1, x2) in resumed.x.iter().zip(&straight.x) {
            assert_eq!(
                x1.to_bits(),
                x2.to_bits(),
                "MINRES bits differ at cut {cut}"
            );
        }
    }
    log("minres-resume", "pass", "3 cut points bitwise");
}

/// Nonsymmetric fixture: convection–diffusion from fs-opdsl on
/// kuhn(2), interior-reduced.
fn convdiff_csr() -> (fs_sparse::Csr, usize) {
    let (complex, positions) = kuhn_cube(2);
    let geo = element_geometry(&complex, &positions);
    let (def, expr) =
        fs_opdsl::fixtures::convection_diffusion(&complex, &geo, 0.5, [2.0, -1.0, 0.5]);
    let full = def.lower(expr).materialize().expect("linear");
    let interior: Vec<usize> = (0..positions.len())
        .filter(|&v| !fs_feec::on_unit_cube_boundary(positions[v]))
        .collect();
    let mut slot = vec![usize::MAX; positions.len()];
    for (i, &v) in interior.iter().enumerate() {
        slot[v] = i;
    }
    let mut red = fs_sparse::Coo::new(interior.len(), interior.len());
    for (i, &v) in interior.iter().enumerate() {
        let (cols, vals) = full.row(v);
        for (&c, &val) in cols.iter().zip(vals) {
            if slot[c] != usize::MAX {
                red.push(i, slot[c], val);
            }
        }
    }
    (red.assemble(), interior.len())
}

#[test]
fn gmres_nonsymmetric_and_transposed() {
    let (a, n) = convdiff_csr();
    let dense = a.to_dense();
    let op = CsrOp::general(a);
    let b = rand_vec(n, 5);
    // Primal.
    let mut st = GmresState::new(&b, 30);
    let rep = st.run(&op, &b, 1e-11, 200, false);
    assert!(rep.converged, "GMRES failed: {rep:?}");
    let lu = fs_la::factor::lu(&dense, n).expect("nonsingular");
    let mut x_ref = b.clone();
    lu.solve(&mut x_ref);
    let dev =
        st.x.iter()
            .zip(&x_ref)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0f64, f64::max);
    let scale = x_ref.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
    assert!(
        dev < 1e-8 * scale.max(1.0),
        "GMRES vs LU deviation {dev:.3e}"
    );
    // TRANSPOSED solve through the SAME machinery (adjoint readiness):
    // Aᵀ y = c, verified against the transposed dense solve.
    let c = rand_vec(n, 6);
    let mut stt = GmresState::new(&c, 30);
    let rept = stt.run(&op, &c, 1e-11, 200, true);
    assert!(rept.converged, "transposed GMRES failed: {rept:?}");
    // Verify: Aᵀ y ≈ c via the operator itself.
    let mut aty = vec![0.0f64; n];
    op.apply_transpose(&stt.x, &mut aty);
    let rel = norm2(
        &c.iter()
            .zip(&aty)
            .map(|(ci, ai)| ci - ai)
            .collect::<Vec<_>>(),
    ) / norm2(&c);
    assert!(rel < 1e-10, "transposed solve residual {rel:.3e}");
    // Comparable convergence (adjoint-readiness gate): within 2x.
    assert!(
        rept.iters <= 2 * rep.iters + 10,
        "transposed convergence degraded: {} vs {}",
        rept.iters,
        rep.iters
    );
    log(
        "gmres",
        "pass",
        &format!("primal iters={} transposed iters={}", rep.iters, rept.iters),
    );
}

#[test]
fn gmres_resume_at_cycle_boundaries() {
    let (a, n) = convdiff_csr();
    let op = CsrOp::general(a);
    let b = rand_vec(n, 7);
    let mut straight = GmresState::new(&b, 10);
    straight.run(&op, &b, 1e-12, 40, false);
    for cut in [1usize, 3] {
        let mut first = GmresState::new(&b, 10);
        first.run(&op, &b, 1e-12, cut, false);
        let mut resumed = first.clone();
        resumed.run(&op, &b, 1e-12, 40 - cut, false);
        assert_eq!(resumed.iters, straight.iters, "iters differ at cut {cut}");
        for (x1, x2) in resumed.x.iter().zip(&straight.x) {
            assert_eq!(x1.to_bits(), x2.to_bits(), "GMRES bits differ at cut {cut}");
        }
    }
    log("gmres-resume", "pass", "cycle-boundary bitwise resume");
}

#[test]
fn stall_diagnosis_is_structured() {
    let (a, n) = poisson_csr(6);
    let op = CsrOp::symmetric(a);
    let b = rand_vec(n, 8);
    // Budget far too small: still falling → BudgetExhausted.
    let mut st = CgState::new(&op, &IdentityPrecond, &b);
    let rep = st.run(&op, &IdentityPrecond, 1e-14, 5);
    assert!(!rep.converged);
    assert_eq!(rep.diagnosis, Some(StallDiagnosis::BudgetExhausted));
    // A genuinely plateauing solve: restarted GMRES on a cyclic
    // shift matrix (the classic stagnation example — with restart
    // m < n every cycle reproduces x = 0 exactly, so the TRUE
    // residual history is flat at 1). CG's recursive residual never
    // plateaus honestly (see diag_probe.rs: singular systems DIVERGE).
    let ns = 24usize;
    let mut shift = fs_sparse::Coo::new(ns, ns);
    for i in 0..ns {
        shift.push(i, (i + 1) % ns, 1.0);
    }
    let op_shift = CsrOp::general(shift.assemble());
    let mut e1 = vec![0.0f64; ns];
    e1[0] = 1.0;
    let mut st2 = GmresState::new(&e1, 5);
    let rep2 = st2.run(&op_shift, &e1, 1e-12, 30, false);
    assert!(!rep2.converged);
    assert_eq!(rep2.diagnosis, Some(StallDiagnosis::Plateau));
    log(
        "diagnosis",
        "pass",
        "BudgetExhausted + Plateau distinguished",
    );
}

#[test]
fn pmg_iteration_counts_flat_across_ladders() {
    // The acceptance gate: CG+pMG iteration counts must stay in a flat
    // envelope across the ORDER ladder and the MESH ladder, while
    // identity-preconditioned counts grow with both.
    let mut table = Vec::new();
    for &(m, r) in &[(2usize, 2usize), (2, 3), (2, 4), (4, 2), (4, 4)] {
        let op = MaskedTensorOp::new(m, r);
        // RANDOM rhs: an eigenfunction load makes ALL Krylov methods
        // converge in a handful of iterations (few effective modes)
        // and would fake both flatness and the identity comparison.
        let mut b = rand_vec(op.n(), 50 + u32::try_from(10 * m + r).expect("small"));
        for (bi, &mk) in b.iter_mut().zip(op.mask()) {
            if !mk {
                *bi = 0.0;
            }
        }
        // Chebyshev degree grows with r (keeps the growth mild).
        let pmg = PMultigrid::new(m, r, r + 2);
        let mut st = CgState::new(&op, &pmg, &b);
        let rep = st.run(&op, &pmg, 1e-10, 100);
        assert!(rep.converged, "pMG-CG failed at m={m} r={r}: {rep:?}");
        let mut st_id = CgState::new(&op, &IdentityPrecond, &b);
        let rep_id = st_id.run(&op, &IdentityPrecond, 1e-10, 5_000);
        table.push((m, r, rep.iters, rep_id.iters));
        log(
            "pmg-ladder",
            "info",
            &format!(
                "m={m} r={r} pmg_iters={} identity_iters={}",
                rep.iters, rep_id.iters
            ),
        );
        assert!(
            rep.iters <= 80,
            "pMG iterations out of envelope at m={m} r={r}: {}",
            rep.iters
        );
        // Solutions agree.
        let dev =
            st.x.iter()
                .zip(&st_id.x)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0f64, f64::max);
        let scale = st_id.x.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
        assert!(
            dev < 1e-7 * scale.max(1.0),
            "pMG solution deviates: {dev:.3e}"
        );
    }
    // "Order-independent-ish", envelope-gated (the bead's words): all
    // pMG counts inside a fixed envelope while identity counts blow up
    // with difficulty — the advantage RATIO at the hard corner is the
    // honest headline number, and mild p-growth is documented rather
    // than hidden (hierarchical-injection transfers; see CONTRACT).
    let pmg_hard = table[4].2;
    let id_hard = table[4].3;
    assert!(
        id_hard >= 10 * pmg_hard,
        "pMG advantage at the hard corner too small: {table:?}"
    );
    let id_easy = table[0].3;
    assert!(
        id_hard > id_easy * 2,
        "identity counts should grow across the ladder: {table:?}"
    );
    log("pmg-envelope", "pass", &format!("{table:?}"));
}

#[test]
fn deterministic_dot_is_length_shaped() {
    // The reduction shape is a function of length only: same values,
    // same result; and it equals the fixed-shape reference combiner.
    let a = rand_vec(1000, 9);
    let b = rand_vec(1000, 10);
    let d1 = dot(&a, &b);
    let d2 = dot(&a, &b);
    assert_eq!(d1.to_bits(), d2.to_bits());
    let prods: Vec<f64> = a.iter().zip(&b).map(|(x, y)| x * y).collect();
    assert_eq!(
        d1.to_bits(),
        fs_tilelang::deterministic_sum(&prods).to_bits()
    );
    log("det-dot", "pass", "fixed-shape reduction");
}

const GOLDEN_HASH: u64 = 0xbc00_5985_1f9c_4a8a; // recorded at tfz.10 slices 1-2, frozen

#[test]
fn solver_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let (a, n) = poisson_csr(4);
    let op = CsrOp::symmetric(a);
    let b = rand_vec(n, 11);
    let mut cg = CgState::new(&op, &IdentityPrecond, &b);
    cg.run(&op, &IdentityPrecond, 1e-11, 2_000);
    for v in cg.x.iter().step_by(3) {
        feed(*v);
    }
    feed(cg.iters as f64);
    let mut mr = MinresState::new(&op, &b);
    mr.run(&op, 1e-10, 2_000);
    for v in mr.x.iter().step_by(5) {
        feed(*v);
    }
    let (c, nc) = convdiff_csr();
    let opc = CsrOp::general(c);
    let bc = rand_vec(nc, 12);
    let mut gm = GmresState::new(&bc, 25);
    gm.run(&opc, &bc, 1e-10, 100, false);
    for v in gm.x.iter().step_by(3) {
        feed(*v);
    }
    // pMG-preconditioned solve output.
    let opt = MaskedTensorOp::new(2, 3);
    let bt = {
        let mut v = rand_vec(opt.n(), 13);
        for (vi, &mk) in v.iter_mut().zip(opt.mask()) {
            if !mk {
                *vi = 0.0;
            }
        }
        v
    };
    let pmg = PMultigrid::new(2, 3, 3);
    let mut st = CgState::new(&opt, &pmg, &bt);
    st.run(&opt, &pmg, 1e-9, 60);
    for v in st.x.iter().step_by(17) {
        feed(*v);
    }
    log("solver-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, GOLDEN_HASH,
        "solver bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}

#[test]
fn mixed_precision_reaches_f64_accuracy() {
    // f32 inner CG + f64 refinement must hit f64-grade residuals on
    // the Poisson fixture (κ ≪ 1/ε_f32), with the work split visible
    // in the report; the total f32 inner iterations replace what a
    // plain f64 CG would spend, at half the memory traffic per
    // iteration (the bandwidth win — measured on the perf lane, not
    // asserted here).
    let (a, n) = poisson_csr(4);
    let a32 = CsrF32::from_csr(&a);
    let op = CsrOp::symmetric(a);
    let b = rand_vec(n, 60);
    let mut x = vec![0.0f64; n];
    let rep = mixed_cg_refine(&op, &a32, &b, &mut x, 1e-12, 20, 2_000);
    assert!(rep.converged, "mixed refinement failed: {rep:?}");
    assert!(rep.rel_residual < 1e-12);
    assert!(
        rep.refine_steps >= 2,
        "expected multiple refinement passes: {rep:?}"
    );
    // f64-grade solution: compare against plain f64 CG.
    let mut st = CgState::new(&op, &IdentityPrecond, &b);
    st.run(&op, &IdentityPrecond, 1e-13, 5_000);
    let dev = x
        .iter()
        .zip(&st.x)
        .map(|(p, q)| (p - q).abs())
        .fold(0.0f64, f64::max);
    let scale = st.x.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
    assert!(
        dev < 1e-9 * scale.max(1.0),
        "mixed vs f64 CG deviation {dev:.3e}"
    );
    log(
        "mixed-precision",
        "pass",
        &format!(
            "refines={} inner_f32_iters={} rel={:.2e}",
            rep.refine_steps, rep.inner_iters, rep.rel_residual
        ),
    );
    // Determinism: repeat run bitwise.
    let mut x2 = vec![0.0f64; n];
    let rep2 = mixed_cg_refine(&op, &a32, &b, &mut x2, 1e-12, 20, 2_000);
    assert_eq!(rep.inner_iters, rep2.inner_iters);
    assert!(x.iter().zip(&x2).all(|(p, q)| p.to_bits() == q.to_bits()));
}

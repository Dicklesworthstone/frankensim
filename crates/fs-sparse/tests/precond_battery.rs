//! Preconditioner battery: spectral-bound enclosure, Chebyshev/ILU(0)
//! acceleration vs plain CG, SA-AMG two-grid quality and
//! grid-independence, anisotropic adversarial fixture, setup determinism,
//! and the cross-ISA golden hash.

use fs_sparse::precond::{
    Chebyshev, IdentityPrecond, Precond, SaAmg, ilu0, lambda_max_estimate, pcg,
};
use fs_sparse::{Coo, Csr};

/// 2D 5-point Laplacian on a g×g grid.
fn laplacian_2d(g: usize) -> Csr {
    let n = g * g;
    let mut coo = Coo::new(n, n);
    for i in 0..g {
        for j in 0..g {
            let u = i * g + j;
            coo.push(u, u, 4.0);
            if i > 0 {
                coo.push(u, u - g, -1.0);
            }
            if i + 1 < g {
                coo.push(u, u + g, -1.0);
            }
            if j > 0 {
                coo.push(u, u - 1, -1.0);
            }
            if j + 1 < g {
                coo.push(u, u + 1, -1.0);
            }
        }
    }
    coo.assemble()
}

/// Anisotropic 2D operator: ε·∂xx + ∂yy (ε ≪ 1 — the classic AMG
/// adversary: smooth error aligns with the strong direction).
fn anisotropic_2d(g: usize, eps: f64) -> Csr {
    let n = g * g;
    let mut coo = Coo::new(n, n);
    for i in 0..g {
        for j in 0..g {
            let u = i * g + j;
            coo.push(u, u, 2.0 * eps + 2.0);
            if j > 0 {
                coo.push(u, u - 1, -eps);
            }
            if j + 1 < g {
                coo.push(u, u + 1, -eps);
            }
            if i > 0 {
                coo.push(u, u - g, -1.0);
            }
            if i + 1 < g {
                coo.push(u, u + g, -1.0);
            }
        }
    }
    coo.assemble()
}

fn rhs(n: usize) -> Vec<f64> {
    // Deterministic, libm-free right-hand side.
    (0..n)
        .map(|i| 1.0 + (((i * 40_503) % 101) as f64) / 101.0 - 0.5)
        .collect()
}

#[test]
fn lambda_max_bound_encloses_analytic_spectrum() {
    // 2D Laplacian spectral max: 4·(1 + cos(π/(g+1))) < 8; the estimate
    // with safety 1.1 must sit at-or-above the true max and below 1.3×.
    let g = 24;
    let a = laplacian_2d(g);
    let true_max = 4.0
        - 2.0 * (std::f64::consts::PI * (g as f64) / ((g + 1) as f64)).cos()
        - 2.0 * (std::f64::consts::PI * (g as f64) / ((g + 1) as f64)).cos();
    let est = lambda_max_estimate(&a, 30, 1.1);
    assert!(
        est >= true_max * 0.999,
        "bound must enclose: est {est} vs true {true_max}"
    );
    assert!(
        est <= true_max * 1.3,
        "bound sloppy: est {est} vs true {true_max}"
    );
    println!(
        "{{\"suite\":\"fs-sparse\",\"case\":\"lambda-max\",\"verdict\":\"pass\",\"detail\":\"est {est:.4} encloses analytic {true_max:.4}\"}}"
    );
}

#[test]
fn chebyshev_and_ilu_accelerate_cg() {
    let g = 32;
    let a = laplacian_2d(g);
    let n = g * g;
    let b = rhs(n);
    let tol = 1e-10;
    // Plain CG baseline.
    let mut x0 = vec![0.0; n];
    let plain = pcg(&a, &b, &mut x0, &IdentityPrecond, tol, 2000);
    assert!(plain.converged, "plain CG must converge: {plain:?}");
    // Chebyshev-preconditioned.
    let cheb = Chebyshev::new(&a, 4, 30.0);
    let mut x1 = vec![0.0; n];
    let c = pcg(&a, &b, &mut x1, &cheb, tol, 2000);
    assert!(c.converged);
    assert!(
        c.iters * 2 < plain.iters,
        "Chebyshev must at least halve iterations: {} vs {}",
        c.iters,
        plain.iters
    );
    // ILU(0)-preconditioned.
    let ilu = ilu0(&a).expect("Laplacian factors without breakdown");
    let mut x2 = vec![0.0; n];
    let i = pcg(&a, &b, &mut x2, &ilu, tol, 2000);
    assert!(i.converged);
    assert!(
        i.iters * 2 < plain.iters,
        "ILU(0) must at least halve iterations: {} vs {}",
        i.iters,
        plain.iters
    );
    // All three agree on the solution.
    for k in 0..n {
        assert!(
            (x1[k] - x0[k]).abs() < 1e-7,
            "chebyshev solution drift at {k}"
        );
        assert!((x2[k] - x0[k]).abs() < 1e-7, "ilu solution drift at {k}");
    }
    println!(
        "{{\"suite\":\"fs-sparse\",\"case\":\"precond-accel\",\"verdict\":\"pass\",\"detail\":\"plain {} vs cheb {} vs ilu {} iters\"}}",
        plain.iters, c.iters, i.iters
    );
}

#[test]
fn ilu_breakdown_is_typed() {
    // Missing diagonal in row 1 → structured breakdown at that row.
    let mut coo = Coo::new(2, 2);
    coo.push(0, 0, 1.0);
    coo.push(0, 1, 1.0);
    coo.push(1, 0, 1.0); // no (1,1) entry
    let a = coo.assemble();
    match ilu0(&a) {
        Err(e) => assert_eq!(e.row, 1, "breakdown row: {e:?}"),
        Ok(_) => panic!("must break down without a diagonal"),
    }
}

#[test]
fn pcg_preserves_rhs_scale_instead_of_squaring_it_away() {
    let a = Csr::identity(2);
    for amplitude in [f64::from_bits(1), 1e-200, 1.0, 1e200, f64::MAX] {
        let b = [amplitude, -amplitude];
        let mut x = [0.0; 2];
        let report = pcg(&a, &b, &mut x, &IdentityPrecond, 1e-12, 2);
        assert!(report.converged, "amplitude {amplitude:e}: {report:?}");
        assert_eq!(x.map(f64::to_bits), b.map(f64::to_bits));
        assert_eq!(report.iters, 1);
        assert_eq!(report.rel_residual, 0.0);
    }
}

#[test]
fn pcg_scaled_curvature_preserves_multistep_solves() {
    let a = Csr::from_parts(3, 3, vec![0, 1, 2, 3], vec![0, 1, 2], vec![2.0, 7.0, 13.0]);
    for amplitude in [1e-200, 1.0, 1e200] {
        let exact = [amplitude, -amplitude, 0.5 * amplitude];
        let mut b = [0.0; 3];
        a.spmv(&exact, &mut b);
        let mut x = [0.0; 3];
        let report = pcg(&a, &b, &mut x, &IdentityPrecond, 1e-12, 8);
        assert!(report.converged, "amplitude {amplitude:e}: {report:?}");
        assert!(
            report.iters > 1,
            "the fixture must exercise the CG recurrence"
        );
        for (&actual, &expected) in x.iter().zip(&exact) {
            assert!((actual / amplitude - expected / amplitude).abs() < 1e-12);
        }
    }
}

#[test]
fn pcg_report_measures_the_returned_iterate() {
    // The recursive FMA residual of this one-step solve is nonzero, whereas
    // independently multiplying the returned x gives exactly the stored RHS.
    let a = Csr::from_parts(1, 1, vec![0, 1], vec![0], vec![3.0]);
    let mut x = [0.0];
    let report = pcg(&a, &[1.0], &mut x, &IdentityPrecond, 0.0, 1);
    let mut ax = [0.0];
    a.spmv(&x, &mut ax);
    assert_eq!(ax, [1.0]);
    assert_eq!(report.rel_residual, (1.0 - ax[0]).abs());
    assert!(report.converged, "{report:?}");
}

#[test]
fn pcg_restarts_when_a_zero_recursive_residual_is_optimistic() {
    // The first step has an exactly zero recursive residual, but rounding
    // 0.1 * 3.0 makes the returned iterate fail the actual equation.
    let a = Csr::from_parts(1, 1, vec![0, 1], vec![0], vec![0.1]);
    let mut capped = [0.0];
    let report = pcg(&a, &[0.3], &mut capped, &IdentityPrecond, 0.0, 1);
    assert!(!report.converged, "{report:?}");
    assert_eq!(report.iters, 1);
    assert!(report.rel_residual > 0.0);
    let mut corrected = [0.0];
    let report = pcg(&a, &[0.3], &mut corrected, &IdentityPrecond, 0.0, 2);
    assert!(report.converged, "{report:?}");
    assert_eq!(report.iters, 2);
    let mut ax = [0.0];
    a.spmv(&corrected, &mut ax);
    assert_eq!(ax, [0.3]);
    assert_eq!(report.rel_residual, 0.0);
}

#[test]
fn pcg_checks_entry_and_iteration_cap_before_preconditioning() {
    struct MustNotApply;
    impl Precond for MustNotApply {
        fn apply(&self, _: &[f64], _: &mut [f64]) {
            panic!("a completed or zero-budget solve must not apply a preconditioner");
        }
    }
    let a = Csr::identity(2);
    let mut exact = [2.0, -3.0];
    let report = pcg(&a, &[2.0, -3.0], &mut exact, &MustNotApply, 1e-12, 10);
    assert!(report.converged);
    assert_eq!(report.iters, 0);
    let mut unchanged = [0.0; 2];
    let report = pcg(&a, &[2.0, -3.0], &mut unchanged, &MustNotApply, 1e-12, 0);
    assert!(!report.converged);
    assert_eq!(report.iters, 0);
    assert_eq!(report.rel_residual, 1.0);
    assert_eq!(unchanged, [0.0; 2]);
}

#[test]
fn pcg_breakdown_keeps_the_last_finite_iterate() {
    for diagonal in [0.0, -1.0] {
        let a = Csr::from_parts(1, 1, vec![0, 1], vec![0], vec![diagonal]);
        let mut x = [0.0];
        let report = pcg(&a, &[1.0], &mut x, &IdentityPrecond, 1e-12, 4);
        assert!(!report.converged, "nonpositive curvature: {report:?}");
        assert_eq!(report.iters, 0);
        assert_eq!(report.rel_residual, 1.0);
        assert_eq!(x, [0.0]);
    }
    struct InvalidPreconditioner;
    impl Precond for InvalidPreconditioner {
        fn apply(&self, _: &[f64], z: &mut [f64]) {
            z.fill(f64::NAN);
        }
    }
    let mut x = [0.0];
    let report = pcg(
        &Csr::identity(1),
        &[1.0],
        &mut x,
        &InvalidPreconditioner,
        1e-12,
        4,
    );
    assert!(!report.converged);
    assert_eq!(report.iters, 0);
    assert_eq!(report.rel_residual, 1.0);
    assert_eq!(x, [0.0]);
}

#[test]
fn pcg_zero_rhs_requires_a_zero_recomputed_residual() {
    let a = Csr::identity(1);
    let mut x = [f64::from_bits(1)];
    let report = pcg(&a, &[0.0], &mut x, &IdentityPrecond, 1e-12, 2);
    assert!(report.converged, "{report:?}");
    assert_eq!(report.iters, 1);
    assert_eq!(report.rel_residual, 0.0);
    assert_eq!(x, [0.0]);
}

#[test]
fn pcg_zero_tolerance_does_not_round_a_nonzero_relative_residual_to_success() {
    let b = [1e300, f64::from_bits(1)];
    let mut x = [1e300, 0.0];
    let report = pcg(&Csr::identity(2), &b, &mut x, &IdentityPrecond, 0.0, 0);
    assert!(
        !report.converged,
        "a finite nonzero residual is not exact convergence"
    );
    assert!(report.rel_residual > 0.0);
    let report = pcg(&Csr::identity(2), &b, &mut x, &IdentityPrecond, 0.0, 1);
    assert!(report.converged, "{report:?}");
    assert_eq!(x.map(f64::to_bits), b.map(f64::to_bits));
}

#[test]
fn amg_grid_independence_and_complexity() {
    let tol = 1e-10;
    let mut iter_counts = Vec::new();
    for g in [32usize, 64] {
        let a = laplacian_2d(g);
        let n = g * g;
        let amg = SaAmg::new(&a, 0.08, 3);
        assert!(
            amg.operator_complexity() < 2.0,
            "operator complexity {} too high at g={g}",
            amg.operator_complexity()
        );
        assert!(
            amg.level_sizes.len() >= 2,
            "hierarchy must coarsen: {:?}",
            amg.level_sizes
        );
        let b = rhs(n);
        let mut x = vec![0.0; n];
        let rep = pcg(&a, &b, &mut x, &amg, tol, 200);
        assert!(rep.converged, "AMG-PCG must converge at g={g}: {rep:?}");
        iter_counts.push(rep.iters);
        println!(
            "{{\"suite\":\"fs-sparse\",\"case\":\"amg\",\"verdict\":\"info\",\"detail\":\"g={g} levels={:?} opcx={:.3} iters={}\"}}",
            amg.level_sizes,
            amg.operator_complexity(),
            rep.iters
        );
    }
    // Grid independence envelope: iteration counts stay low and close.
    assert!(
        iter_counts.iter().all(|&c| c <= 40),
        "AMG iteration counts blew up: {iter_counts:?}"
    );
    let diff = iter_counts[1].abs_diff(iter_counts[0]);
    assert!(
        diff <= 10,
        "iterations should be near grid-independent: {iter_counts:?}"
    );
    println!(
        "{{\"suite\":\"fs-sparse\",\"case\":\"amg-independence\",\"verdict\":\"pass\",\"detail\":\"iters {iter_counts:?}\"}}"
    );
}

#[test]
fn amg_survives_anisotropy() {
    // ε = 1e-3: strength-based aggregation must follow the strong (y)
    // direction; convergence may degrade but must not fail.
    let g = 32;
    let a = anisotropic_2d(g, 1e-3);
    let n = g * g;
    let amg = SaAmg::new(&a, 0.08, 3);
    let b = rhs(n);
    let mut x = vec![0.0; n];
    let rep = pcg(&a, &b, &mut x, &amg, 1e-9, 300);
    assert!(
        rep.converged,
        "anisotropic fixture must still converge: {rep:?}"
    );
    println!(
        "{{\"suite\":\"fs-sparse\",\"case\":\"amg-anisotropy\",\"verdict\":\"pass\",\"detail\":\"eps=1e-3 iters={} levels={:?}\"}}",
        rep.iters, amg.level_sizes
    );
}

#[test]
fn setup_and_solve_are_deterministic() {
    let g = 24;
    let a = laplacian_2d(g);
    let n = g * g;
    let run = || {
        let amg = SaAmg::new(&a, 0.08, 3);
        let b = rhs(n);
        let mut x = vec![0.0; n];
        let rep = pcg(&a, &b, &mut x, &amg, 1e-10, 200);
        (amg.level_sizes.clone(), rep.iters, x)
    };
    let (l1, i1, x1) = run();
    let (l2, i2, x2) = run();
    assert_eq!(l1, l2, "hierarchy must be rerun-identical");
    assert_eq!(i1, i2, "iteration count must be rerun-identical");
    assert!(
        x1.iter().zip(&x2).all(|(p, q)| p.to_bits() == q.to_bits()),
        "solutions must be BITWISE identical"
    );
}

/// Recorded on aarch64-apple (M4 Pro); must match on x86-64 (trj).
const GOLDEN_HASH: u64 = 0x752f_215a_26e3_2fea;

#[test]
fn precond_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let g = 24;
    let a = laplacian_2d(g);
    let n = g * g;
    let b = rhs(n);
    // Chebyshev apply bits.
    let cheb = Chebyshev::new(&a, 4, 30.0);
    let mut z = vec![0.0; n];
    cheb.apply(&b, &mut z);
    for &v in &z {
        feed(v);
    }
    // ILU(0)-PCG solution bits.
    let ilu = ilu0(&a).unwrap();
    let mut x = vec![0.0; n];
    pcg(&a, &b, &mut x, &ilu, 1e-10, 500);
    for &v in &x {
        feed(v);
    }
    // AMG-PCG solution bits + hierarchy shape.
    let amg = SaAmg::new(&a, 0.08, 3);
    let mut y = vec![0.0; n];
    let rep = pcg(&a, &b, &mut y, &amg, 1e-10, 200);
    for &v in &y {
        feed(v);
    }
    feed(rep.iters as f64);
    for &s in &amg.level_sizes {
        feed(s as f64);
    }
    println!(
        "{{\"suite\":\"fs-sparse\",\"case\":\"precond-golden\",\"verdict\":\"info\",\"detail\":\"{acc:#018x}\"}}"
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "preconditioner bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with \
         semantic justification (golden-evidence policy)"
    );
}

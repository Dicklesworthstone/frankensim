//! E2E casebook for the sparse direct facility (bead
//! frankensim-fsim-sparse-direct-4a38j): factor a stiffness pencil
//! (K − σM) at a ladder of five shifts with ONE symbolic analysis, emit
//! JSON-line evidence per factorization, and certify every inertia triple
//! against the ANALYTIC eigenvalue counts of the fixture.
//!
//! Fixture: K = 5-point Laplacian on an s×s grid (the plate-pencil stand-in
//! until the plates bead lands), M = identity mass. Eigenvalues are known in
//! closed form: λ_{ij} = (2 − 2cos(iπ/(s+1))) + (2 − 2cos(jπ/(s+1))), so the
//! inertia of K − σM is checkable without any numeric eigensolver — the
//! spectrum-slicing certification pattern the vibration-eig bead consumes.

use fs_sparse::{Coo, Csr, DirectOrdering, LdltOptions, SymbolicLdlt};

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

/// K − σ·I with the SAME pattern as K (the diagonal is already stored).
fn shifted(k: &Csr, sigma: f64) -> Csr {
    let n = k.nrows();
    let mut coo = Coo::new(n, n);
    for r in 0..n {
        let (cols, vals) = k.row(r);
        for (&c, &v) in cols.iter().zip(vals) {
            coo.push(r, c, v);
        }
    }
    for i in 0..n {
        coo.push(i, i, -sigma);
    }
    coo.assemble()
}

/// Analytic count of grid-Laplacian eigenvalues strictly below `sigma`.
fn analytic_count_below(s: usize, sigma: f64) -> usize {
    let mut count = 0;
    for i in 1..=s {
        for j in 1..=s {
            let li = 2.0 - 2.0 * (i as f64 * std::f64::consts::PI / (s as f64 + 1.0)).cos();
            let lj = 2.0 - 2.0 * (j as f64 * std::f64::consts::PI / (s as f64 + 1.0)).cos();
            if li + lj < sigma {
                count += 1;
            }
        }
    }
    count
}

#[test]
fn shift_ladder_casebook_with_certified_inertia() {
    let s = 100; // 10_000 unknowns
    let n = s * s;
    let k = laplacian_2d(s);
    let t0 = std::time::Instant::now();
    let sym = SymbolicLdlt::analyze(&k, DirectOrdering::Amd).expect("analyze");
    let analyze_ms = t0.elapsed().as_secs_f64() * 1e3;
    // Fill comparison against the natural-ordering baseline (recorded, and
    // bounded: AMD must not be worse — the scrambled-ordering guard).
    let nat = SymbolicLdlt::analyze(&k, DirectOrdering::Natural).expect("analyze natural");
    assert!(
        sym.nnz_l() <= nat.nnz_l(),
        "AMD fill {} exceeded natural-ordering fill {}",
        sym.nnz_l(),
        nat.nnz_l()
    );
    println!(
        "{{\"suite\":\"fs-sparse-direct-casebook\",\"stage\":\"analyze\",\"n\":{n},\"nnz_a\":{},\"nnz_l_amd\":{},\"nnz_l_natural\":{},\"supernodes\":{},\"analyze_ms\":{analyze_ms:.1}}}",
        k.nnz(),
        sym.nnz_l(),
        nat.nnz_l(),
        sym.supernode_count()
    );

    // Interior shifts spanning the spectrum (0, 8); K − σI is INDEFINITE at
    // every interior shift — this exercises exactly what shift-invert needs.
    let shifts = [0.05, 0.5, 2.0, 4.0, 7.0];
    let mut prev_below = 0usize;
    for &sigma in &shifts {
        let a = shifted(&k, sigma);
        let t1 = std::time::Instant::now();
        let f = sym.factor(&a, &LdltOptions::default()).expect("factor");
        let factor_ms = t1.elapsed().as_secs_f64() * 1e3;
        let inertia = f.inertia();
        let below = analytic_count_below(s, sigma);
        assert_eq!(
            inertia.negative, below,
            "inertia(negative) must equal the analytic eigenvalue count below sigma={sigma}"
        );
        assert_eq!(inertia.positive, n - below, "inertia must sum to n");
        assert!(
            inertia.negative >= prev_below,
            "eigenvalue counts must be monotone in sigma"
        );
        prev_below = inertia.negative;

        // Solve gate at each shift: manufactured solution, residual bound.
        let x_true: Vec<f64> = (0..n).map(|i| ((i % 17) as f64) - 8.0).collect();
        let mut b = vec![0.0; n];
        a.spmv(&x_true, &mut b);
        let x = f.solve(&b);
        let mut ax = vec![0.0; n];
        a.spmv(&x, &mut ax);
        let xin = x.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
        let resid = ax
            .iter()
            .zip(&b)
            .fold(0.0f64, |m, (&p, &q)| m.max((p - q).abs()));
        let tol = 1024.0 * (n as f64) * f64::EPSILON * 8.0 * xin.max(1.0);
        assert!(
            resid <= tol,
            "residual {resid} above gate {tol} at sigma={sigma}"
        );

        let st = f.stats();
        println!(
            "{{\"suite\":\"fs-sparse-direct-casebook\",\"stage\":\"factor\",\"sigma\":{sigma},\"nnz_a\":{},\"nnz_l\":{},\"fill_ratio\":{:.3},\"flops\":{},\"peak_front_bytes\":{},\"max_front_dim\":{},\"pivots_1x1\":{},\"pivots_2x2\":{},\"inertia_pos\":{},\"inertia_neg\":{},\"analytic_below\":{below},\"solve_resid\":{resid:.3e},\"factor_ms\":{factor_ms:.1}}}",
            st.nnz_a,
            st.nnz_l,
            st.fill_ratio,
            st.flops,
            st.peak_front_bytes,
            st.max_front_dim,
            st.pivots_1x1,
            st.pivots_2x2,
            inertia.positive,
            inertia.negative,
        );
    }
}

#[test]
fn repeat_factorization_is_bitwise_identical_at_scale() {
    // G5-style determinism audit on the casebook fixture (smaller grid so
    // the double factorization stays cheap in debug CI).
    let k = laplacian_2d(24);
    let a = shifted(&k, 1.5);
    let sym = SymbolicLdlt::analyze(&k, DirectOrdering::Amd).expect("analyze");
    let f1 = sym.factor(&a, &LdltOptions::default()).expect("factor 1");
    let f2 = sym.factor(&a, &LdltOptions::default()).expect("factor 2");
    let b: Vec<f64> = (0..k.nrows()).map(|i| ((i % 13) as f64) - 6.0).collect();
    let x1 = f1.solve(&b);
    let x2 = f2.solve(&b);
    assert!(
        x1.iter().zip(&x2).all(|(p, q)| p.to_bits() == q.to_bits()),
        "repeat factorization+solve must be bitwise identical"
    );
    println!(
        "{{\"suite\":\"fs-sparse-direct-casebook\",\"stage\":\"determinism\",\"n\":{},\"verdict\":\"bitwise\"}}",
        k.nrows()
    );
}

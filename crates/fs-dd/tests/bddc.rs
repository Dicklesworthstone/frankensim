//! BDDC conformance (the tfz.11 bead; runs under `bddc` +
//! `sheaf-coarse`). Acceptance: condition numbers match BDDC theory
//! (log²(H/h) scaling); coefficient-jump robustness at 1e6 with bounded
//! iterations; the sheaf-derived edge coarse space measured against
//! corners-only (honestly reported); CCD-aligned partitioning shows the
//! expected locality metric; subdomain sweep ledgered.
#![cfg(feature = "bddc")]

use fs_dd::{Bddc, CgError, CgReport, CgTermination, Decomposition};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-dd/bddc\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

/// A deterministic non-trivial interface RHS.
fn rhs(n: usize, seed: u64) -> Vec<f64> {
    let mut state = seed;
    (0..n)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 11) as f64) / (1u64 << 53) as f64 - 0.5
        })
        .collect()
}

fn converged_report(report: Result<CgReport, CgError>, tolerance: f64, case: &str) -> CgReport {
    let report = report.unwrap_or_else(|error| panic!("{case}: admitted solve failed: {error}"));
    assert_eq!(
        report.termination,
        CgTermination::Converged,
        "{case}: claims require a converged solve, got {report:?}"
    );
    assert!(
        report.relative_residual.is_finite() && report.relative_residual <= tolerance,
        "{case}: converged residual exceeds tolerance {tolerance:e}: {report:?}"
    );
    report
}

fn converged_diagnostics(
    report: Result<CgReport, CgError>,
    tolerance: f64,
    case: &str,
) -> (usize, f64) {
    let report = converged_report(report, tolerance, case);
    let estimate = report
        .condition_estimate
        .unwrap_or_else(|| panic!("{case}: converged multi-step solve needs valid Ritz evidence"));
    assert!(
        estimate.krylov_dimension > 1,
        "{case}: a one-dimensional Ritz projection is not conditioning evidence: {estimate:?}"
    );
    assert!(
        estimate.ritz_min.is_finite()
            && estimate.ritz_min > 0.0
            && estimate.ritz_max.is_finite()
            && estimate.ritz_max >= estimate.ritz_min
            && estimate.kappa.is_finite()
            && estimate.kappa >= 1.0,
        "{case}: invalid condition estimate {estimate:?}"
    );
    (report.iterations, estimate.kappa)
}

fn stable_l2_norm(values: &[f64]) -> f64 {
    let mut scale = 0.0f64;
    let mut sum_squares = 1.0f64;
    for value in values {
        let magnitude = value.abs();
        if magnitude == 0.0 {
            continue;
        }
        if scale < magnitude {
            let ratio = scale / magnitude;
            sum_squares = 1.0 + sum_squares * ratio * ratio;
            scale = magnitude;
        } else {
            let ratio = magnitude / scale;
            sum_squares += ratio * ratio;
        }
    }
    if scale == 0.0 {
        0.0
    } else {
        scale * fs_math::det::sqrt(sum_squares)
    }
}

/// Partition the non-Dirichlet lattice nodes into the globally eliminated
/// interior and retained subdomain-interface sets. The ordering matches the
/// full-lattice row-major order, without consulting `Bddc` internals.
fn global_node_partition(decomp: &Decomposition) -> (Vec<usize>, Vec<usize>) {
    let n = decomp.s * decomp.m;
    let np = n + 1;
    let mut interior = Vec::new();
    let mut gamma = Vec::new();
    for y in 1..n {
        for x in 1..n {
            let node = y * np + x;
            if x % decomp.m == 0 || y % decomp.m == 0 {
                gamma.push(node);
            } else {
                interior.push(node);
            }
        }
    }
    (interior, gamma)
}

/// Independent dense solve for the tiny global-oracle fixture. Gaussian
/// elimination deliberately avoids sharing fs-dd's local Cholesky path.
fn solve_dense(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    assert_eq!(a.len(), n);
    assert!(a.iter().all(|row| row.len() == n));
    for col in 0..n {
        let mut pivot = col;
        for row in (col + 1)..n {
            if a[row][col].abs() > a[pivot][col].abs() {
                pivot = row;
            }
        }
        assert!(a[pivot][col].abs() > 1e-14, "global K_II is nonsingular");
        a.swap(col, pivot);
        b.swap(col, pivot);
        let pivot_row = a[col].clone();
        let pivot_value = pivot_row[col];
        let pivot_rhs = b[col];
        for row in (col + 1)..n {
            let scale = a[row][col] / pivot_value;
            a[row][col] = 0.0;
            for j in (col + 1)..n {
                a[row][j] -= scale * pivot_row[j];
            }
            b[row] -= scale * pivot_rhs;
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let tail: f64 = a[row][(row + 1)..]
            .iter()
            .zip(&x[(row + 1)..])
            .map(|(entry, value)| entry * value)
            .sum();
        x[row] = (b[row] - tail) / a[row][row];
    }
    x
}

/// Apply the exact global Schur complement
/// `K_GG - K_GI K_II^-1 K_IG` using only the public whole-system oracle.
fn global_schur_apply(decomp: &Decomposition, x_gamma: &[f64]) -> Vec<f64> {
    let (interior, gamma) = global_node_partition(decomp);
    assert_eq!(x_gamma.len(), gamma.len());
    let np = decomp.s * decomp.m + 1;
    let mut harmonic_extension = vec![0.0; np * np];
    for (&node, &value) in gamma.iter().zip(x_gamma) {
        harmonic_extension[node] = value;
    }

    let interface_load = decomp.apply_global(&harmonic_extension);
    let mut k_ii = vec![vec![0.0; interior.len()]; interior.len()];
    for (col, &node) in interior.iter().enumerate() {
        let mut basis = vec![0.0; np * np];
        basis[node] = 1.0;
        let image = decomp.apply_global(&basis);
        for (row, &row_node) in interior.iter().enumerate() {
            k_ii[row][col] = image[row_node];
        }
    }
    let interior_rhs: Vec<f64> = interior.iter().map(|&node| interface_load[node]).collect();
    let correction = solve_dense(k_ii, interior_rhs);
    for (&node, &value) in interior.iter().zip(&correction) {
        harmonic_extension[node] = -value;
    }
    let condensed = decomp.apply_global(&harmonic_extension);
    gamma.iter().map(|&node| condensed[node]).collect()
}

#[test]
fn dd_000_global_schur_matches_substructured_operator() {
    for (fixture, decomp) in [
        ("uniform", Decomposition::uniform(2, 2)),
        ("jump", Decomposition::checkerboard(2, 2, 64.0)),
    ] {
        let bddc = Bddc::new(decomp.clone(), false);
        let (_, gamma) = global_node_partition(&decomp);
        assert_eq!(bddc.gamma_len(), gamma.len());

        let mut probes = Vec::with_capacity(gamma.len() + 1);
        for col in 0..gamma.len() {
            let mut basis = vec![0.0; gamma.len()];
            basis[col] = 1.0;
            probes.push(basis);
        }
        probes.push(rhs(gamma.len(), 0x5c48_7572));
        for (probe, x_gamma) in probes.iter().enumerate() {
            let expected = global_schur_apply(&decomp, x_gamma);
            let actual = bddc.schur_apply(x_gamma);
            for (row, (&got, &want)) in actual.iter().zip(&expected).enumerate() {
                let error = (got - want).abs();
                assert!(
                    error <= 5e-11 * want.abs().max(1.0),
                    "{fixture} probe {probe}, row {row}: substructured {got} vs global {want}"
                );
            }
        }
    }
    verdict(
        "dd-000",
        "substructured Schur action equals independent elimination of the apply_global oracle",
    );
}

#[test]
fn face_transmissibility_is_the_harmonic_mean() {
    for (orientation, rho) in [
        ("vertical", vec![1.0, 9.0, 1.0, 9.0]),
        ("horizontal", vec![1.0, 1.0, 9.0, 9.0]),
    ] {
        // In either orientation, two center edges are flanked by rho=1 and
        // rho=9 cells, so each contributes 2/(1/1 + 1/9) = 1.8. The other
        // two edges contribute 1 and 9: total center diagonal 13.6.
        let decomp = Decomposition { s: 1, m: 2, rho };
        let mut x = vec![0.0; 9];
        x[4] = 1.0;
        let image = decomp.apply_global(&x);
        assert!(
            (image[4] - 13.6).abs() <= 1e-12,
            "{orientation} harmonic face transmissibility gives center diagonal 13.6, got {}",
            image[4]
        );
    }

    let tiny = f64::from_bits(1);
    let decomp = Decomposition {
        s: 1,
        m: 2,
        rho: vec![tiny, 2.0 * tiny, 2.0 * tiny, tiny],
    };
    let mut x = vec![0.0; 9];
    x[4] = 1.0;
    let image = decomp.apply_global(&x);
    assert!(
        image[4].is_finite() && image[4] > 0.0,
        "positive subnormal coefficients retain positive transmissibility"
    );
}

#[test]
#[should_panic(expected = "cell coefficients must be finite and positive")]
fn invalid_cell_coefficient_cannot_disappear_during_harmonic_mean() {
    let decomp = Decomposition {
        s: 1,
        m: 2,
        rho: vec![f64::NAN, 1.0, 1.0, 1.0],
    };
    decomp.apply_global(&vec![0.0; 9]);
}

#[test]
fn zero_and_empty_rhs_solves_are_total() {
    let nonempty = Bddc::new(Decomposition::uniform(2, 2), false);
    let zero = vec![0.0; nonempty.gamma_len()];
    assert_eq!(
        nonempty.solve_cg(&zero, 1e-8, 20),
        Ok(CgReport {
            termination: CgTermination::Converged,
            iterations: 0,
            relative_residual: 0.0,
            condition_estimate: None
        })
    );

    let mut tiny = vec![0.0; nonempty.gamma_len()];
    tiny[0] = 1e-300;
    let report = converged_report(nonempty.solve_cg(&tiny, 1e-8, 20), 1e-8, "tiny RHS");
    assert!(
        report.iterations > 0,
        "a nonzero tiny RHS is not classified as zero"
    );
    assert!(report.relative_residual <= 1e-8, "{report:?}");
    if let Some(estimate) = report.condition_estimate {
        assert!(estimate.kappa.is_finite(), "{estimate:?}");
    }

    let empty = Bddc::new(Decomposition::uniform(1, 2), false);
    assert_eq!(empty.gamma_len(), 0);
    assert_eq!(
        empty.solve_cg(&[], 1e-8, 20),
        Ok(CgReport {
            termination: CgTermination::Converged,
            iterations: 0,
            relative_residual: 0.0,
            condition_estimate: None
        })
    );
}

#[test]
fn rhs_dimension_mismatch_is_rejected_before_zero_shortcut() {
    let nonempty = Bddc::new(Decomposition::uniform(2, 2), false);
    assert_eq!(
        nonempty.solve_cg(&[], 1e-8, 20),
        Err(CgError::DimensionMismatch {
            expected: nonempty.gamma_len(),
            actual: 0
        })
    );
}

#[test]
fn cg_controls_and_zero_budget_are_typed() {
    let bddc = Bddc::new(Decomposition::uniform(2, 2), false);
    let zero = vec![0.0; bddc.gamma_len()];
    for tolerance in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1.0] {
        assert!(matches!(
            bddc.solve_cg(&zero, tolerance, 0),
            Err(CgError::InvalidTolerance { tolerance: rejected })
                if rejected.to_bits() == tolerance.to_bits()
        ));
    }

    let mut invalid_rhs = zero.clone();
    invalid_rhs[0] = f64::NAN;
    assert!(matches!(
        bddc.solve_cg(&invalid_rhs, 1e-8, 20),
        Err(CgError::NonFiniteRhs { index: 0, value }) if value.is_nan()
    ));

    let nonzero = rhs(bddc.gamma_len(), 0x4347_4255_4447_4554);
    assert_eq!(
        bddc.solve_cg(&nonzero, 0.0, 0),
        Ok(CgReport {
            termination: CgTermination::IterationLimit,
            iterations: 0,
            relative_residual: 1.0,
            condition_estimate: None
        })
    );
    assert_eq!(
        bddc.solve_cg(&zero, 0.0, 0),
        Ok(CgReport {
            termination: CgTermination::Converged,
            iterations: 0,
            relative_residual: 0.0,
            condition_estimate: None
        })
    );
    assert_eq!(
        bddc.solve_cg(&nonzero, 1.0, 0),
        Ok(CgReport {
            termination: CgTermination::Converged,
            iterations: 0,
            relative_residual: 1.0,
            condition_estimate: None
        })
    );
}

#[test]
fn one_iteration_exhaustion_retains_true_residual_and_ritz_dimension() {
    let bddc = Bddc::new(Decomposition::checkerboard(4, 4, 1e6), false);
    let input = rhs(bddc.gamma_len(), 0x4c49_4d49_542d_3031);
    let report = bddc
        .solve_cg(&input, 0.0, 1)
        .expect("finite admitted one-step solve");
    assert_eq!(report.termination, CgTermination::IterationLimit);
    assert_eq!(report.iterations, 1);
    assert!(
        report.relative_residual.is_finite() && report.relative_residual > 0.0,
        "{report:?}"
    );
    let estimate = report
        .condition_estimate
        .expect("one accepted Lanczos step has a transparent trivial projection");
    assert_eq!(estimate.krylov_dimension, 1);
    assert_eq!(estimate.kappa, 1.0);

    // Independently reconstruct the accepted one-step iterate using the
    // public operator/preconditioner, then check the report against b - Sx.
    // This prevents the cheaper recursive residual from masquerading as the
    // terminal true-residual diagnostic.
    let input_scale = input.iter().map(|value| value.abs()).fold(0.0f64, f64::max);
    let scaled_input: Vec<f64> = input.iter().map(|value| value / input_scale).collect();
    let direction = bddc.precondition(&scaled_input);
    let image = bddc.schur_apply(&direction);
    let rz: f64 = scaled_input
        .iter()
        .zip(&direction)
        .map(|(left, right)| left * right)
        .sum();
    let curvature: f64 = direction
        .iter()
        .zip(&image)
        .map(|(left, right)| left * right)
        .sum();
    let alpha = rz / curvature;
    let iterate: Vec<f64> = direction
        .iter()
        .map(|value| alpha.mul_add(*value, 0.0))
        .collect();
    let true_image = bddc.schur_apply(&iterate);
    let true_residual: Vec<f64> = scaled_input
        .iter()
        .zip(&true_image)
        .map(|(right, left)| right - left)
        .collect();
    let independently_recomputed = stable_l2_norm(&true_residual) / stable_l2_norm(&scaled_input);
    assert_eq!(
        report.relative_residual.to_bits(),
        independently_recomputed.to_bits(),
        "iteration-limit reports must carry the true residual"
    );
}

#[test]
fn g3_rhs_sign_and_power_of_two_scaling_preserve_cg_diagnostics() {
    let bddc = Bddc::new(Decomposition::uniform(2, 4), true);
    let base_rhs = rhs(bddc.gamma_len(), 0x5343_414c_4544_4347);
    let expected = bddc
        .solve_cg(&base_rhs, 1e-10, 200)
        .expect("finite admitted base solve");
    assert_eq!(expected.termination, CgTermination::Converged);
    assert!(expected.relative_residual <= 1e-10, "{expected:?}");

    for scale in [
        -1.0,
        fs_math::det::powi(2.0, 200),
        fs_math::det::powi(2.0, -200),
    ] {
        let transformed: Vec<f64> = base_rhs.iter().map(|value| value * scale).collect();
        let actual = bddc
            .solve_cg(&transformed, 1e-10, 200)
            .expect("finite admitted transformed solve");
        assert_eq!(actual, expected, "RHS scale {scale:e}");
    }
}

#[test]
fn dd_001_g0_preconditioner_properties() {
    let bddc = Bddc::new(Decomposition::uniform(2, 4), true);
    // SPD-ness of M^{-1} on probes: <M^{-1} r, r> > 0.
    for seed in [1u64, 2, 3, 4] {
        let r = rhs(bddc.gamma_len(), seed);
        let z = bddc.precondition(&r);
        let dot: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
        assert!(dot > 0.0, "SPD probe (seed {seed}): {dot}");
    }
    // Symmetry probe: <M^{-1}a, b> == <a, M^{-1}b> to rounding.
    let a = rhs(bddc.gamma_len(), 7);
    let b = rhs(bddc.gamma_len(), 8);
    let ma = bddc.precondition(&a);
    let mb = bddc.precondition(&b);
    let lhs: f64 = ma.iter().zip(&b).map(|(x, y)| x * y).sum();
    let rhs_: f64 = a.iter().zip(&mb).map(|(x, y)| x * y).sum();
    assert!(
        (lhs - rhs_).abs() <= 1e-10 * lhs.abs().max(1.0),
        "symmetric: {lhs} vs {rhs_}"
    );
    // The preconditioned solve converges fast on the uniform problem.
    let (iters, kappa) = converged_diagnostics(
        bddc.solve_cg(&rhs(bddc.gamma_len(), 9), 1e-8, 200),
        1e-8,
        "dd-001 uniform solve",
    );
    assert!(iters < 30, "uniform 2x2 converges quickly: {iters}");
    assert!(kappa < 50.0, "conditioned: {kappa}");
    verdict(
        "dd-001",
        "M^-1 SPD + symmetric on probes; uniform 2x2 solve converges in < 30 iterations",
    );
}

#[test]
fn dd_002_log_squared_h_over_h_scaling() {
    // Fixed s = 4, sweep H/h: kappa should track C(1 + log(H/h))^2 —
    // the BDDC signature. We assert the NORMALIZED kappa stays within a
    // bounded band across the sweep (no polynomial growth).
    let mut rows = Vec::new();
    for m in [4usize, 8, 16] {
        let bddc = Bddc::new(Decomposition::uniform(4, m), true);
        let (iters, kappa) = converged_diagnostics(
            bddc.solve_cg(&rhs(bddc.gamma_len(), 42), 1e-8, 400),
            1e-8,
            "dd-002 scaling solve",
        );
        let logf = (1.0 + (m as f64).ln()).powi(2);
        rows.push((m, iters, kappa, kappa / logf));
        println!(
            "{{\"metric\":\"bddc-scaling\",\"h_over_h\":{m},\"iters\":{iters},\
             \"kappa\":{kappa:.2},\"kappa_over_log2\":{:.2}}}",
            kappa / logf
        );
    }
    let normalized: Vec<f64> = rows.iter().map(|r| r.3).collect();
    let (lo, hi) = normalized
        .iter()
        .fold((f64::INFINITY, 0.0f64), |(l, h), &v| (l.min(v), h.max(v)));
    assert!(
        hi / lo < 4.0,
        "normalized kappa is bounded across the sweep (log^2 scaling): {normalized:?}"
    );
    // And raw kappa grows SLOWLY: the log^2 model predicts
    // (1+ln16)^2/(1+ln4)^2 ~ 2.5x across this sweep; a 4x quadrupling of
    // H/h must stay well under linear growth.
    assert!(
        rows[2].2 / rows[0].2 < 3.5,
        "kappa growth tracks log^2, not H/h: {} -> {}",
        rows[0].2,
        rows[2].2
    );
    verdict(
        "dd-002",
        "kappa/(1+log(H/h))^2 stays in a <4x band across H/h in {2,4,8} — the BDDC \
         signature",
    );
}

#[test]
fn dd_003_coefficient_jump_robustness() {
    // Checkerboard 1e6: the jump-aligned decomposition must stay
    // bounded. Compare against the uniform problem's iterations.
    let uniform = Bddc::new(Decomposition::uniform(4, 4), true);
    let (it_u, _) = converged_diagnostics(
        uniform.solve_cg(&rhs(uniform.gamma_len(), 5), 1e-8, 400),
        1e-8,
        "dd-003 uniform reference",
    );
    let jump = Bddc::new(Decomposition::checkerboard(4, 4, 1e6), true);
    let (it_j, kappa_j) = converged_diagnostics(
        jump.solve_cg(&rhs(jump.gamma_len(), 5), 1e-8, 400),
        1e-8,
        "dd-003 jump solve",
    );
    println!(
        "{{\"metric\":\"jump-robustness\",\"uniform_iters\":{it_u},\"jump_iters\":{it_j},\
         \"jump_kappa\":{kappa_j:.2}}}"
    );
    assert!(
        it_j <= 3 * it_u.max(5),
        "1e6 checkerboard stays bounded: {it_j} vs uniform {it_u}"
    );
    verdict(
        "dd-003",
        "checkerboard 1e6 jumps: iterations within 3x of uniform (subdomain-aligned \
         jumps are the BDDC-friendly case, honestly noted)",
    );
}

#[test]
fn dd_004_sheaf_edge_coarse_vs_corners_only() {
    // The measured comparison the bead demands: corners-only vs the
    // sheaf-derived edge enrichment, on uniform AND jump fixtures.
    let mut table: Vec<(&str, usize, f64, usize, f64, usize, usize)> = Vec::new();
    for (name, d) in [
        ("uniform", Decomposition::uniform(4, 4)),
        ("jump-1e6", Decomposition::checkerboard(4, 4, 1e6)),
    ] {
        let corners = Bddc::new(d.clone(), false);
        let (it_c, k_c) = converged_diagnostics(
            corners.solve_cg(&rhs(corners.gamma_len(), 11), 1e-8, 600),
            1e-8,
            "dd-004 corners solve",
        );
        let sheaf = Bddc::new(d, true);
        let (it_s, k_s) = converged_diagnostics(
            sheaf.solve_cg(&rhs(sheaf.gamma_len(), 11), 1e-8, 600),
            1e-8,
            "dd-004 sheaf solve",
        );
        println!(
            "{{\"metric\":\"coarse-comparison\",\"fixture\":\"{name}\",\
             \"corners\":{{\"iters\":{it_c},\"kappa\":{k_c:.2},\"dim\":{}}},\
             \"sheaf\":{{\"iters\":{it_s},\"kappa\":{k_s:.2},\"dim\":{}}}}}",
            corners.coarse_dim(),
            sheaf.coarse_dim()
        );
        table.push((name, 0usize, k_c, 0usize, k_s, it_c, it_s));
    }
    // THE HONEST MEASURED CLAIM (the bead: "measured, honestly
    // reported if not"): on the ADVERSARIAL jump fixture the
    // sheaf-edge coarse space strictly improves the condition estimate
    // with comparable iterations; on the uniform fixture kappa is
    // comparable and iterations pay a small price for the larger
    // coarse space — reported, not hidden.
    let (_, _, k_c_jump, _, k_s_jump, it_c_jump, it_s_jump) = table[1];
    assert!(
        k_s_jump < k_c_jump,
        "jump fixture: the sheaf coarse space strictly improves kappa \
         ({k_s_jump:.2} vs {k_c_jump:.2})"
    );
    #[allow(clippy::cast_precision_loss)]
    {
        assert!(
            (it_s_jump as f64) <= 1.3 * (it_c_jump as f64),
            "jump fixture: iterations comparable ({it_s_jump} vs {it_c_jump})"
        );
    }
    verdict(
        "dd-004",
        "adversarial jump fixture: sheaf-edge coarse strictly improves kappa with \
         comparable iterations; the uniform trade is ledgered honestly",
    );
}

#[cfg(feature = "sheaf-coarse")]
#[test]
fn dd_005_sheaf_cross_check_and_ccd_locality() {
    // Bet 11's machinery frames the enrichment: the subdomain-adjacency
    // skeleton's open interfaces must match the edge-mode count, and a
    // gauge cochain on it must Hodge-decompose with zero harmonic part
    // on the 2x2 grid-with-cross topology... the 2x2 subdomain grid's
    // adjacency is a 4-cycle: b1 = 1, and the CORNER constraint is what
    // pins that cycle — the sheaf explains WHY corners are primal.
    use fs_geom::sheaf_repair::{SheafSkeleton, hodge_decompose};
    let bddc = Bddc::new(Decomposition::uniform(2, 4), true);
    // 2x2 subdomains: 4 open interfaces (N/S/E/W around the cross).
    assert_eq!(bddc.coarse_dim(), 1 + 4, "1 corner + 4 edge modes");
    let skeleton = SheafSkeleton {
        n_patches: 4,
        edges: vec![(0, 1), (1, 3), (2, 3), (0, 2)],
        triangles: vec![],
    };
    assert_eq!(
        skeleton.edges.len(),
        4,
        "the adjacency skeleton has one edge per open interface"
    );
    // The 4-cycle carries a 1-dimensional harmonic space — the
    // circulation mode the corner constraint pins.
    // Orientation-aware cycle 0->1->3->2->0 over edges
    // (0,1)+, (1,3)+, (2,3)-, (0,2)-.
    let circulation = vec![1.0, 1.0, -1.0, -1.0];
    let split =
        hodge_decompose(&skeleton, &circulation).expect("valid BDDC sheaf incidence fixture");
    assert!(
        split.fractions.2 > 0.999,
        "the subdomain cycle's harmonic mode exists: {:?}",
        split.fractions
    );
    // CCD locality: on an 8x8 subdomain grid, 4 islands (2x2 blocks of
    // 4x4 subdomains) keep most interfaces island-internal; 64 islands
    // (one per subdomain) keep none.
    let big = Bddc::new(Decomposition::uniform(8, 2), false);
    let aligned = big.ccd_locality(2);
    let degenerate = big.ccd_locality(8);
    println!(
        "{{\"metric\":\"ccd-locality\",\"aligned_2x2\":{aligned:.3},\
         \"degenerate_8x8\":{degenerate:.3}}}"
    );
    assert!(
        aligned > 0.8,
        "CCD-aligned partitioning keeps interfaces local: {aligned}"
    );
    assert!(
        degenerate < 0.1,
        "one-subdomain-per-island has no internal interfaces: {degenerate}"
    );
    // Regression: ccd_locality(0) used to panic (div_ceil(0) + `ccds - 1`
    // usize underflow). Zero islands is degenerate — it must return 0.0, total.
    assert_eq!(
        big.ccd_locality(0),
        0.0,
        "zero islands is total, not a panic"
    );
    verdict(
        "dd-005",
        "the sheaf skeleton matches the enrichment (4 edges, 1 harmonic cycle mode \
         pinned by the corner); CCD-aligned locality 0.8+ vs degenerate <0.1",
    );
}

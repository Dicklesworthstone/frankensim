//! Resumable eigensolver-service battery (bead bfid). Printed
//! measurements on every gate; the dense reference FALSIFIES the
//! sparse-path claims, it never proves them.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_spectral::service::{
    CertifiedEigenvalue, DenseSymOp, EigenBackend, EigenQuery, EigenService, ServiceError,
    SymmetricOp, gap_report,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xBF1D,
                kernel_id: 13,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

/// Deterministic seeded symmetric matrix with a controlled spectrum:
/// Q^T diag(spectrum) Q for a fixed Givens-rotation cascade Q.
fn seeded_sym(spectrum: &[f64]) -> DenseSymOp {
    let n = spectrum.len();
    let mut a = vec![0.0f64; n * n];
    for (i, &s) in spectrum.iter().enumerate() {
        a[i * n + i] = s;
    }
    // Apply a fixed cascade of exact Givens rotations G a Gᵀ so the
    // matrix is dense but the spectrum is EXACTLY the input list
    // (rotations are similarity transforms).
    let mut rotate = |p: usize, q: usize, c: f64, s: f64, a: &mut Vec<f64>| {
        for j in 0..n {
            let (apj, aqj) = (a[p * n + j], a[q * n + j]);
            a[p * n + j] = c * apj - s * aqj;
            a[q * n + j] = s * apj + c * aqj;
        }
        for i in 0..n {
            let (aip, aiq) = (a[i * n + p], a[i * n + q]);
            a[i * n + p] = c * aip - s * aiq;
            a[i * n + q] = s * aip + c * aiq;
        }
    };
    let (c, s) = (0.6, 0.8); // exact 3-4-5 rotation
    for k in 0..(n - 1) {
        rotate(k, (k + 3) % n, c, s, &mut a);
    }
    // Symmetrize exactly against accumulated roundoff so the operator
    // constructor's exact-symmetry gate passes.
    for i in 0..n {
        for j in (i + 1)..n {
            let m = 0.5 * (a[i * n + j] + a[j * n + i]);
            a[i * n + j] = m;
            a[j * n + i] = m;
        }
    }
    DenseSymOp::new(n, a).expect("seeded operator is square, finite, symmetric")
}

fn spectrum_a() -> Vec<f64> {
    (0..48).map(|i| -3.0 + 0.25 * i as f64).collect()
}

#[test]
fn sv_001_accuracy_vs_dense_reference_both_backends_both_ends() {
    with_cx(|cx| {
        let spectrum = spectrum_a();
        let op = seeded_sym(&spectrum);
        let mut sorted = spectrum.clone();
        sorted.sort_by(f64::total_cmp);
        for backend in [EigenBackend::Lanczos, EigenBackend::Lobpcg] {
            for largest in [true, false] {
                let mut service = EigenService::new(
                    backend,
                    op.dim(),
                    EigenQuery {
                        k: 4,
                        largest,
                        tol: 1e-9,
                        steps_per_tick: 8,
                    },
                )
                .expect("query is valid");
                let progress = service
                    .run_to_tolerance(&op, cx, 64)
                    .expect("converges on a 48-dof dense fixture");
                let truth: Vec<f64> = if largest {
                    sorted.iter().rev().take(4).copied().collect()
                } else {
                    sorted.iter().take(4).copied().collect()
                };
                for t in &truth {
                    let hit = progress
                        .pairs
                        .iter()
                        .any(|p| p.interval.0 - 1e-8 <= *t && *t <= p.interval.1 + 1e-8);
                    assert!(
                        hit,
                        "{backend:?} largest={largest}: true eigenvalue {t} missed by every \
                         certified interval"
                    );
                }
                println!(
                    "sv-001 {backend:?} largest={largest}: {} pairs in {} ticks, worst residual {:e}",
                    progress.pairs.len(),
                    progress.ticks,
                    progress
                        .pairs
                        .iter()
                        .map(|p| p.residual)
                        .fold(0.0f64, f64::max)
                );
            }
        }
    });
}

#[test]
fn sv_002_degenerate_spectrum_clusters_at_resolution() {
    with_cx(|cx| {
        // Fourfold-degenerate top eigenvalue, well separated from the
        // rest.
        let mut spectrum = vec![2.0, 2.0, 2.0, 2.0];
        spectrum.extend((0..44).map(|i| -4.0 + 0.1 * i as f64));
        let op = seeded_sym(&spectrum);
        let mut service = EigenService::new(
            EigenBackend::Lobpcg,
            op.dim(),
            EigenQuery {
                k: 4,
                largest: true,
                tol: 1e-8,
                steps_per_tick: 8,
            },
        )
        .expect("query is valid");
        let progress = service
            .run_to_tolerance(&op, cx, 96)
            .expect("degenerate block converges");
        let report = gap_report(&progress.pairs);
        println!(
            "sv-002: clusters {:?} leading gap lower bound {:e}",
            report
                .clusters
                .iter()
                .map(|c| (c.hull, c.count))
                .collect::<Vec<_>>(),
            report.leading_gap_lower_bound
        );
        // All four pairs cluster at the degenerate value at this
        // resolution ("count at resolution", not exact multiplicity).
        let top = report
            .clusters
            .iter()
            .find(|c| c.hull.0 <= 2.0 && 2.0 <= c.hull.1)
            .expect("a cluster covers the degenerate eigenvalue");
        assert_eq!(top.count, 4, "fourfold cluster at resolution");
    });
}

#[test]
fn sv_003_split_run_is_bitwise_equal_to_straight_run() {
    with_cx(|cx| {
        let op = seeded_sym(&spectrum_a());
        let query = EigenQuery {
            k: 3,
            largest: true,
            tol: 1e-30, // never converges: exercise raw ticking
            steps_per_tick: 4,
        };
        for backend in [EigenBackend::Lanczos, EigenBackend::Lobpcg] {
            let mut straight = EigenService::new(backend, op.dim(), query).expect("valid");
            for _ in 0..6 {
                let _ = straight.tick(&op, cx).expect("tick");
            }
            let final_straight = straight.tick(&op, cx).expect("tick");

            let mut split = EigenService::new(backend, op.dim(), query).expect("valid");
            for _ in 0..3 {
                let _ = split.tick(&op, cx).expect("tick");
            }
            // Checkpoint = plain clone; resume from the checkpoint.
            let mut resumed = split.clone();
            drop(split);
            for _ in 0..3 {
                let _ = resumed.tick(&op, cx).expect("tick");
            }
            let final_split = resumed.tick(&op, cx).expect("tick");

            assert_eq!(final_straight.pairs.len(), final_split.pairs.len());
            for (a, b) in final_straight.pairs.iter().zip(final_split.pairs.iter()) {
                assert_eq!(
                    a.value.to_bits(),
                    b.value.to_bits(),
                    "{backend:?}: split-run eigenvalue bits differ"
                );
                assert_eq!(
                    a.residual.to_bits(),
                    b.residual.to_bits(),
                    "{backend:?}: split-run residual bits differ"
                );
            }
            println!("sv-003 {backend:?}: split == straight over 7 ticks (bitwise)");
        }
    });
}

#[test]
fn sv_004_warm_start_speeds_up_continuation() {
    with_cx(|cx| {
        // A(θ) = A0 + θ·D with D a small seeded symmetric perturbation.
        let base = spectrum_a();
        let op0 = seeded_sym(&base);
        let perturbed: Vec<f64> = base.iter().map(|s| s * 1.02 + 0.01).collect();
        let op1 = seeded_sym(&perturbed);
        let query = EigenQuery {
            k: 3,
            largest: true,
            tol: 1e-9,
            steps_per_tick: 2,
        };
        let mut cold0 = EigenService::new(EigenBackend::Lobpcg, op0.dim(), query).expect("valid");
        let sol0 = cold0
            .run_to_tolerance(&op0, cx, 200)
            .expect("theta0 converges");
        // Continuation: warm-start theta1 from theta0's Ritz vectors.
        let seeds: Vec<Vec<f64>> = sol0
            .pairs
            .iter()
            .rev() // largest first, matching the block's target end
            .take(3)
            .map(|p| p.vector.clone())
            .collect();
        let mut warm = EigenService::warm(EigenBackend::Lobpcg, op1.dim(), query, &seeds)
            .expect("warm seed accepted");
        let warm_sol = warm
            .run_to_tolerance(&op1, cx, 200)
            .expect("warm converges");
        let mut cold1 = EigenService::new(EigenBackend::Lobpcg, op1.dim(), query).expect("valid");
        let cold_sol = cold1
            .run_to_tolerance(&op1, cx, 200)
            .expect("cold converges");
        println!(
            "sv-004: warm {} ticks vs cold {} ticks (logged perf evidence, no absolute claim)",
            warm_sol.ticks, cold_sol.ticks
        );
        assert!(
            warm_sol.ticks <= cold_sol.ticks,
            "warm start must not be slower on a small continuation step: warm {} cold {}",
            warm_sol.ticks,
            cold_sol.ticks
        );
    });
}

#[test]
fn sv_005_typed_refusals() {
    with_cx(|cx| {
        let op = seeded_sym(&spectrum_a());
        let bad_k = EigenService::new(
            EigenBackend::Lanczos,
            op.dim(),
            EigenQuery {
                k: 0,
                largest: true,
                tol: 1e-9,
                steps_per_tick: 4,
            },
        );
        assert!(matches!(bad_k, Err(ServiceError::InvalidQuery { .. })));
        let bad_tol = EigenService::new(
            EigenBackend::Lanczos,
            op.dim(),
            EigenQuery {
                k: 2,
                largest: true,
                tol: 0.0,
                steps_per_tick: 4,
            },
        );
        assert!(matches!(bad_tol, Err(ServiceError::InvalidQuery { .. })));
        let bad_block = EigenService::new(
            EigenBackend::Lobpcg,
            6,
            EigenQuery {
                k: 3,
                largest: true,
                tol: 1e-9,
                steps_per_tick: 4,
            },
        );
        assert!(matches!(bad_block, Err(ServiceError::InvalidQuery { .. })));
        let bad_steps = EigenService::new(
            EigenBackend::Lanczos,
            op.dim(),
            EigenQuery {
                k: 2,
                largest: true,
                tol: 1e-9,
                steps_per_tick: 0,
            },
        );
        assert!(matches!(bad_steps, Err(ServiceError::InvalidQuery { .. })));
        let bad_seed = EigenService::warm(
            EigenBackend::Lanczos,
            op.dim(),
            EigenQuery {
                k: 2,
                largest: true,
                tol: 1e-9,
                steps_per_tick: 4,
            },
            &[vec![f64::NAN; op.dim()]],
        );
        assert!(matches!(bad_seed, Err(ServiceError::InvalidSeed)));
        // Rank-deficient warm block: two identical seed vectors.
        let same = vec![1.0f64; op.dim()];
        let rank_deficient = EigenService::warm(
            EigenBackend::Lobpcg,
            op.dim(),
            EigenQuery {
                k: 2,
                largest: true,
                tol: 1e-9,
                steps_per_tick: 4,
            },
            &[same.clone(), same],
        );
        assert!(
            matches!(rank_deficient, Err(ServiceError::InvalidSeed)),
            "rank-deficient warm block must refuse, not be silently completed"
        );
        // Dimension mismatch between operator and service.
        let mut svc = EigenService::new(
            EigenBackend::Lanczos,
            op.dim() + 1,
            EigenQuery {
                k: 2,
                largest: true,
                tol: 1e-9,
                steps_per_tick: 4,
            },
        )
        .expect("valid");
        assert!(matches!(
            svc.tick(&op, cx),
            Err(ServiceError::DimensionMismatch { .. })
        ));
        // Budget exhaustion is typed and resumable.
        let mut tiny_budget = EigenService::new(
            EigenBackend::Lanczos,
            op.dim(),
            EigenQuery {
                k: 3,
                largest: true,
                tol: 1e-30,
                steps_per_tick: 1,
            },
        )
        .expect("valid");
        let out = tiny_budget.run_to_tolerance(&op, cx, 2);
        assert!(
            matches!(out, Err(ServiceError::Unconverged { ticks: 2, .. })),
            "budget exhaustion must be a typed, resumable outcome"
        );
        assert!(tiny_budget.tick(&op, cx).is_ok());
        println!("sv-005: every refusal typed; service resumable after budget exhaustion");
    });
}

#[test]
fn sv_006_intervals_contain_reference_eigenvalues() {
    with_cx(|cx| {
        // Falsification: dense reference values must land inside the
        // certified Weyl intervals of converged pairs.
        let spectrum = spectrum_a();
        let op = seeded_sym(&spectrum);
        let mut service = EigenService::new(
            EigenBackend::Lanczos,
            op.dim(),
            EigenQuery {
                k: 5,
                largest: false,
                tol: 1e-10,
                steps_per_tick: 8,
            },
        )
        .expect("valid");
        let progress = service.run_to_tolerance(&op, cx, 64).expect("converges");
        let mut sorted = spectrum;
        sorted.sort_by(f64::total_cmp);
        for (pair, truth) in progress.pairs.iter().zip(sorted.iter()) {
            println!(
                "sv-006: interval [{:.12}, {:.12}] vs reference {truth:.12}",
                pair.interval.0, pair.interval.1
            );
            assert!(
                pair.interval.0 - 1e-9 <= *truth && *truth <= pair.interval.1 + 1e-9,
                "reference eigenvalue {truth} escapes the certified interval"
            );
        }
    });
}

#[test]
fn sv_007_hostile_sizes_refuse_without_overflow() {
    with_cx(|cx| {
        // fs-la seeded-start constructors: checked arithmetic on
        // public usize inputs — hostile sizes refuse, never wrap or
        // panic.
        use fs_la::eigen::LobpcgState;
        assert!(LobpcgState::with_block(usize::MAX, 2, &[]).is_none());
        assert!(LobpcgState::with_block(usize::MAX, usize::MAX / 2, &[]).is_none());
        assert!(LobpcgState::with_block(12, 0, &[]).is_none());
        assert!(
            LobpcgState::with_block(12, 5, &[0.0; 60]).is_none(),
            "3b > n refuses"
        );
        // Service-level warm seed with hostile dimensions refuses
        // typed.
        let seeds = vec![vec![0.0f64; 8]];
        let out = EigenService::warm(
            EigenBackend::Lobpcg,
            usize::MAX,
            EigenQuery {
                k: 2,
                largest: true,
                tol: 1e-9,
                steps_per_tick: 1,
            },
            &seeds,
        );
        assert!(out.is_err(), "hostile warm dimensions must refuse");
        let _ = cx;
        println!("sv-007: hostile sizes refuse without overflow or allocation");
    });
}

#[test]
fn sv_008_identity_operator_exhausts_subspace_typed() {
    with_cx(|cx| {
        // The identity's Krylov space from any start vector is
        // one-dimensional: asking Lanczos for k = 2 must surface a
        // TYPED subspace-exhausted outcome — never a zero-vector
        // re-entry, division by zero, or NaN pair (the persistent-
        // breakdown regression).
        let n = 16;
        let mut a = vec![0.0f64; n * n];
        for i in 0..n {
            a[i * n + i] = 1.0;
        }
        let op = DenseSymOp::new(n, a).expect("identity is symmetric");
        let mut service = EigenService::new(
            EigenBackend::Lanczos,
            n,
            EigenQuery {
                k: 2,
                largest: true,
                tol: 1e-12,
                steps_per_tick: 4,
            },
        )
        .expect("valid");
        let out = service.run_to_tolerance(&op, cx, 16);
        match out {
            Err(ServiceError::SubspaceExhausted { available }) => {
                println!("sv-008: typed subspace exhaustion with {available} pair(s)");
                assert_eq!(available, 1, "identity yields exactly one Krylov direction");
            }
            other => panic!("sv-008: expected SubspaceExhausted, got {other:?}"),
        }
        // Resumable and stable: further ticks keep returning finite
        // pairs from the exhausted subspace instead of corrupting it.
        let progress = service.tick(&op, cx).expect("post-exhaustion tick is safe");
        assert!(progress.subspace_exhausted);
        assert!(progress.pairs.iter().all(|p| p.value.is_finite()));
        assert!((progress.pairs[0].value - 1.0).abs() < 1e-10);
    });
}

#[test]
fn sv_009_gap_report_merges_bridging_intervals() {
    // Regression for the clustering defect: [0,0], [1,1] look like two
    // clusters until the bridging interval [-1,5] arrives; sorting by
    // interval LOWER bound and extending BOTH hull endpoints must
    // collapse all three into ONE cluster with no false leading gap.
    let mk = |value: f64, lo: f64, hi: f64| CertifiedEigenvalue {
        value,
        residual: (hi - lo) / 2.0,
        interval: (lo, hi),
        vector: Vec::new(),
    };
    let pairs = vec![mk(0.0, 0.0, 0.0), mk(1.0, 1.0, 1.0), mk(2.0, -1.0, 5.0)];
    let report = gap_report(&pairs);
    println!(
        "sv-009: clusters {:?} leading gap {:e}",
        report
            .clusters
            .iter()
            .map(|c| (c.hull, c.count))
            .collect::<Vec<_>>(),
        report.leading_gap_lower_bound
    );
    assert_eq!(
        report.clusters.len(),
        1,
        "bridging interval merges all three"
    );
    assert_eq!(report.clusters[0].count, 3);
    assert_eq!(report.clusters[0].hull, (-1.0, 5.0));
    assert_eq!(report.leading_gap_lower_bound, 0.0);
    // Non-finite intervals never reach clustering.
    let with_nan = vec![mk(0.0, f64::NAN, f64::NAN), mk(1.0, 1.0, 1.0)];
    let filtered = gap_report(&with_nan);
    assert_eq!(filtered.clusters.len(), 1);
    assert_eq!(filtered.clusters[0].count, 1);
}

#[test]
fn sv_010_sparse_operator_seam_via_fs_sparse() {
    with_cx(|cx| {
        // The bead's sparse seam: a CSR-backed operator adapts to
        // SymmetricOp at L1 (fs-sparse is L1) and the service resolves
        // its extremal eigenvalues, falsified against the dense
        // reference on the same matrix.
        struct CsrSymOp {
            n: usize,
            csr: fs_sparse::Csr,
        }
        impl SymmetricOp for CsrSymOp {
            fn dim(&self) -> usize {
                self.n
            }
            fn apply(&self, x: &[f64], y: &mut [f64]) {
                self.csr.spmv(x, y);
            }
        }
        // Symmetric tridiagonal Laplacian: eigenvalues are known
        // analytically: 2 - 2 cos(j*pi/(n+1)).
        let n = 40;
        let mut coo = fs_sparse::Coo::new(n, n);
        for i in 0..n {
            coo.push(i, i, 2.0);
            if i + 1 < n {
                coo.push(i, i + 1, -1.0);
                coo.push(i + 1, i, -1.0);
            }
        }
        let op = CsrSymOp {
            n,
            csr: coo.assemble(),
        };
        let mut service = EigenService::new(
            EigenBackend::Lanczos,
            n,
            EigenQuery {
                k: 3,
                largest: true,
                tol: 1e-9,
                steps_per_tick: 8,
            },
        )
        .expect("valid");
        let progress = service.run_to_tolerance(&op, cx, 64).expect("converges");
        for (rank, pair) in progress.pairs.iter().rev().enumerate() {
            let j = n - rank;
            let truth =
                2.0 - 2.0 * fs_math::det::cos(j as f64 * std::f64::consts::PI / (n as f64 + 1.0));
            println!(
                "sv-010: sparse pair {} interval [{:.12}, {:.12}] vs analytic {truth:.12}",
                rank, pair.interval.0, pair.interval.1
            );
            assert!(
                pair.interval.0 - 1e-8 <= truth && truth <= pair.interval.1 + 1e-8,
                "analytic eigenvalue {truth} escapes the certified interval"
            );
        }
    });
}

//! fs-cheb budget/admission battery (bead frankensim-sj31i.55, slice 1).
//!
//! G0 boundary tables for every admission cap and checked size
//! formula (including `usize::MAX`-shaped requests that must refuse
//! BEFORE allocating), typed refusals where the classic APIs panic,
//! bitwise parity between the budgeted and classic paths, real
//! cancellation with deterministic RESUME equivalence, and receipt
//! determinism (G5).

use asupersync::types::Budget;
use fs_cheb::{
    BuildRun, Cheb1, ChebBudget, ChebError, EigsRun, admit_adaptive_build, admit_dirichlet_eigs,
    admit_root_scan, dirichlet_laplace_eigs_budgeted, try_build_budgeted,
};
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};

fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    if cancelled {
        gate.request();
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x55,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

/// cb-001 — admission boundaries: each cap refuses one-over and admits
/// at-cap; `usize::MAX`-shaped requests refuse via CHECKED formulas
/// before any allocation (no panic, no OOM, no saturation loop).
#[test]
fn cb_001_admission_boundaries() {
    // Adaptive build: coefficients cap exactly at the final grid.
    let mut b = ChebBudget::default();
    b.max_coefficients = 1024;
    admit_adaptive_build(0.0, 1.0, 1024, 16, &b).expect("at-cap admits");
    assert!(matches!(
        admit_adaptive_build(0.0, 1.0, 1025, 16, &b),
        Err(ChebError::CapExceeded {
            what: "retained coefficients",
            ..
        })
    ));
    // Samples cap: total = 2 * final grid.
    let mut b = ChebBudget::default();
    b.max_samples = 2047;
    assert!(matches!(
        admit_adaptive_build(0.0, 1.0, 1024, 16, &b),
        Err(ChebError::CapExceeded {
            what: "adaptive samples",
            ..
        })
    ));
    b.max_samples = 2048;
    admit_adaptive_build(0.0, 1.0, 1024, 16, &b).expect("exact sample budget admits");
    // Temp bytes cap: 2 * grid * 8.
    let mut b = ChebBudget::default();
    b.max_temp_bytes = 2 * 1024 * 8 - 1;
    assert!(matches!(
        admit_adaptive_build(0.0, 1.0, 1024, 16, &b),
        Err(ChebError::CapExceeded {
            what: "adaptive temporary bytes",
            ..
        })
    ));
    // usize::MAX degree: checked next_power_of_two refuses, no alloc.
    assert!(matches!(
        admit_adaptive_build(0.0, 1.0, usize::MAX, 16, &ChebBudget::default()),
        Err(ChebError::Overflow { .. })
    ));
    // Eigensolve: usize::MAX dimension refuses before allocation.
    assert!(matches!(
        admit_dirichlet_eigs(usize::MAX, 1, &ChebBudget::default()),
        Err(ChebError::Overflow { .. })
    ));
    // Eigensolve dimension cap boundary (m = n + 1).
    let mut b = ChebBudget::default();
    b.max_eigen_dim = 64;
    admit_dirichlet_eigs(63, 3, &b).expect("m = 64 admits");
    assert!(matches!(
        admit_dirichlet_eigs(64, 3, &b),
        Err(ChebError::CapExceeded {
            what: "collocation dimension",
            ..
        })
    ));
    // Shape refusals: degenerate eigensolves are caller bugs, not work.
    assert!(matches!(
        admit_dirichlet_eigs(1, 1, &ChebBudget::default()),
        Err(ChebError::Shape { .. })
    ));
    assert!(matches!(
        admit_dirichlet_eigs(24, 0, &ChebBudget::default()),
        Err(ChebError::Shape { .. })
    ));
    assert!(matches!(
        admit_dirichlet_eigs(24, 24, &ChebBudget::default()),
        Err(ChebError::Shape { .. })
    ));
    // Root scan: usize::MAX-coefficient scan refuses via checked mul.
    assert!(matches!(
        admit_root_scan(usize::MAX, &ChebBudget::default()),
        Err(ChebError::Overflow { .. })
    ));
    // Zero caps refuse everything (deterministic first violation).
    let mut zero = ChebBudget::default();
    zero.max_coefficients = 0;
    assert!(admit_adaptive_build(0.0, 1.0, 16, 16, &zero).is_err());
}

/// cb-002 — domain refusals are typed, not panics: NaN, infinite, and
/// reversed endpoints all refuse with the endpoint bits named.
#[test]
fn cb_002_domain_refusals() {
    for (a, b) in [
        (f64::NAN, 1.0),
        (0.0, f64::INFINITY),
        (1.0, 1.0),
        (2.0, 1.0),
    ] {
        let refusal = admit_adaptive_build(a, b, 64, 16, &ChebBudget::default())
            .expect_err("invalid domain must refuse");
        assert!(matches!(refusal, ChebError::Domain { .. }), "{refusal}");
    }
}

/// cb-003 — bitwise parity: on the happy path the budgeted entry
/// points produce EXACTLY the classic results (same sample sequence,
/// same transforms, same plateau/truncation, same eigenvalue bits).
#[test]
fn cb_003_budgeted_matches_classic_bitwise() {
    with_cx(false, |cx| {
        let f = |x: f64| fs_math::det::sin(3.0 * x) + 0.25 * fs_math::det::cos(11.0 * x);
        let classic = Cheb1::build(&f, -1.0, 2.0, 4096);
        let run = try_build_budgeted(&f, -1.0, 2.0, 4096, None, &ChebBudget::default(), cx)
            .expect("budgeted build admits");
        let BuildRun::Complete { function, receipt } = run else {
            panic!("uncancelled build must complete");
        };
        assert_eq!(function.domain(), classic.domain());
        assert_eq!(function.coeffs().len(), classic.coeffs().len());
        for (lhs, rhs) in function.coeffs().iter().zip(classic.coeffs()) {
            assert_eq!(lhs.to_bits(), rhs.to_bits(), "coefficient bit parity");
        }
        assert!(receipt.rounds_completed >= 1 && receipt.samples_spent >= 16);

        let classic_eigs = fs_cheb::dirichlet_laplace_eigs(24, 3);
        let run = dirichlet_laplace_eigs_budgeted(24, 3, &ChebBudget::default(), cx)
            .expect("budgeted eigensolve admits");
        let EigsRun::Complete { eigs, .. } = run else {
            panic!("uncancelled eigensolve must complete");
        };
        assert_eq!(eigs.len(), classic_eigs.len());
        for (lhs, rhs) in eigs.iter().zip(&classic_eigs) {
            assert_eq!(lhs.to_bits(), rhs.to_bits(), "eigenvalue bit parity");
        }

        let poly = Cheb1::build(&|x: f64| (x - 0.25) * (x + 0.5), -1.0, 1.0, 64);
        let classic_roots = poly.roots();
        let budgeted_roots = poly
            .roots_budgeted(&ChebBudget::default(), cx)
            .expect("root scan admits");
        assert_eq!(budgeted_roots.len(), classic_roots.len());
        for (lhs, rhs) in budgeted_roots.iter().zip(&classic_roots) {
            assert_eq!(lhs.to_bits(), rhs.to_bits(), "root bit parity");
        }
    });
}

/// cb-004 — cancellation and RESUME: a pre-cancelled gate drains at
/// the first bounded boundary with an explicit Cancelled state (and a
/// resume point for the constructor); resuming completes with results
/// bitwise-identical to the uncancelled run.
#[test]
fn cb_004_cancellation_and_resume() {
    let f = |x: f64| fs_math::det::exp(-x * x) * fs_math::det::sin(9.0 * x);
    let cancelled = with_cx(true, |cx| {
        try_build_budgeted(&f, -1.0, 1.0, 4096, None, &ChebBudget::default(), cx)
            .expect("admission precedes cancellation")
    });
    let BuildRun::Cancelled {
        resume_from,
        receipt,
    } = cancelled
    else {
        panic!("pre-cancelled gate must drain, not complete");
    };
    assert_eq!(resume_from, 16, "drains before the first round");
    assert_eq!(receipt.samples_spent, 0, "no work after the drain point");

    let resumed = with_cx(false, |cx| {
        try_build_budgeted(
            &f,
            -1.0,
            1.0,
            4096,
            Some(resume_from),
            &ChebBudget::default(),
            cx,
        )
        .expect("resume admits")
    });
    let direct = with_cx(false, |cx| {
        try_build_budgeted(&f, -1.0, 1.0, 4096, None, &ChebBudget::default(), cx)
            .expect("direct admits")
    });
    let (BuildRun::Complete { function: a, .. }, BuildRun::Complete { function: b, .. }) =
        (resumed, direct)
    else {
        panic!("both runs complete");
    };
    for (lhs, rhs) in a.coeffs().iter().zip(b.coeffs()) {
        assert_eq!(lhs.to_bits(), rhs.to_bits(), "resume is bitwise-equivalent");
    }

    // Eigensolve: pre-cancelled drains with an EMPTY converged prefix.
    let run = with_cx(true, |cx| {
        dirichlet_laplace_eigs_budgeted(24, 3, &ChebBudget::default(), cx)
            .expect("admission precedes cancellation")
    });
    let EigsRun::Cancelled { partial_eigs, .. } = run else {
        panic!("pre-cancelled eigensolve must drain");
    };
    assert!(partial_eigs.is_empty(), "no shift completed");

    // Root scan: cancellation refuses with NO partial claim.
    let poly = Cheb1::build(&|x: f64| (x - 0.25) * (x + 0.5), -1.0, 1.0, 64);
    let refusal = with_cx(true, |cx| {
        poly.roots_budgeted(&ChebBudget::default(), cx)
            .expect_err("cancelled scan refuses")
    });
    assert!(matches!(refusal, ChebError::Cancelled), "{refusal}");
}

/// cb-005 — typed refusals where the classic API panics: an
/// unresolvable (discontinuous) function and a non-finite sample both
/// come back as errors from the budgeted path.
#[test]
fn cb_005_typed_refusals_replace_panics() {
    with_cx(false, |cx| {
        let step = |x: f64| if x < 0.5 { -1.0 } else { 1.0 };
        let refusal = try_build_budgeted(&step, 0.0, 1.0, 128, None, &ChebBudget::default(), cx)
            .expect_err("a step function cannot reach the plateau");
        assert!(
            matches!(refusal, ChebError::Unresolved { max_degree: 128 }),
            "{refusal}"
        );

        let singular = |x: f64| 1.0 / (x - 0.5);
        let refusal =
            try_build_budgeted(&singular, 0.0, 1.0, 128, None, &ChebBudget::default(), cx)
                .expect_err("a pole inside the domain cannot sample finitely");
        assert!(
            matches!(
                refusal,
                ChebError::NonFinite { .. } | ChebError::Unresolved { .. }
            ),
            "{refusal}"
        );
    });
}

/// cb-006 — receipt determinism (G5): identical budgeted runs produce
/// identical receipts and identical terminal states.
#[test]
fn cb_006_receipt_determinism() {
    let f = |x: f64| fs_math::det::cos(5.0 * x);
    let run = || {
        with_cx(false, |cx| {
            try_build_budgeted(&f, 0.0, 3.0, 2048, None, &ChebBudget::default(), cx)
                .expect("admits")
        })
    };
    assert_eq!(run(), run(), "whole-run determinism incl. receipts");

    let eig_run = || {
        with_cx(false, |cx| {
            dirichlet_laplace_eigs_budgeted(16, 2, &ChebBudget::default(), cx).expect("admits")
        })
    };
    assert_eq!(
        eig_run(),
        eig_run(),
        "eigensolve determinism incl. receipts"
    );
}

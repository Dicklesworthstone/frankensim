//! fs-cheb budget/admission battery (bead frankensim-sj31i.55, slice 1).
//!
//! G0 boundary tables for every admission cap and checked size
//! formula (including `usize::MAX`-shaped requests that must refuse
//! BEFORE allocating), typed refusals where the classic APIs panic,
//! bitwise parity between the budgeted and classic paths, real
//! cancellation with deterministic RESUME equivalence, and receipt
//! deterministic replay on one admitted execution profile.

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

fn with_gate_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x55,
                kernel_id: 8,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx)
    })
}

/// cb-001 — admission boundaries: each cap refuses one-over and admits
/// at-cap; `usize::MAX`-shaped requests refuse via CHECKED formulas
/// before any allocation (no panic, no OOM, no saturation loop).
#[test]
#[allow(clippy::too_many_lines)] // one ordered boundary table preserves refusal precedence
fn cb_001_admission_boundaries() {
    // Adaptive build: coefficients cap exactly at the final grid.
    let mut b = ChebBudget::default();
    b.max_coefficients = 1024;
    admit_adaptive_build(0.0, 1.0, 1024, 16, &b).expect("at-cap admits");
    assert!(matches!(
        admit_adaptive_build(0.0, 1.0, 2048, 16, &b),
        Err(ChebError::CapExceeded {
            what: "retained coefficients",
            ..
        })
    ));
    // Samples cap: exact geometric sum = 2 * final grid - start.
    let mut b = ChebBudget::default();
    b.max_samples = 2031;
    assert!(matches!(
        admit_adaptive_build(0.0, 1.0, 1024, 16, &b),
        Err(ChebError::CapExceeded {
            what: "adaptive samples",
            ..
        })
    ));
    b.max_samples = 2032;
    admit_adaptive_build(0.0, 1.0, 1024, 16, &b).expect("exact sample budget admits");
    // Temp bytes cap includes sampling plus the complete DCT-II envelope.
    let mut b = ChebBudget::default();
    b.max_temp_bytes = 64 * 1024 - 1;
    assert!(matches!(
        admit_adaptive_build(0.0, 1.0, 1024, 16, &b),
        Err(ChebError::CapExceeded {
            what: "adaptive temporary bytes",
            ..
        })
    ));
    // usize::MAX degree: the conservative cap or checked arithmetic refuses,
    // depending on the caller's explicit cap schedule; neither path allocates.
    assert!(matches!(
        admit_adaptive_build(0.0, 1.0, usize::MAX, 16, &ChebBudget::default()),
        Err(ChebError::CapExceeded { .. }) | Err(ChebError::Overflow { .. })
    ));
    // Eigensolve: usize::MAX dimension refuses before allocation.
    assert!(matches!(
        admit_dirichlet_eigs(usize::MAX, 1, &ChebBudget::default()),
        Err(ChebError::Overflow { .. })
    ));
    let mut enormous = ChebBudget::default();
    enormous.max_eigen_dim = usize::MAX;
    enormous.max_temp_bytes = u64::MAX;
    enormous.max_work_ops = u64::MAX;
    assert!(matches!(
        admit_dirichlet_eigs(usize::MAX - 1, 1, &enormous),
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
    // Conservative eigensolve work/temp envelopes admit exactly at their
    // reported bound and reject one below it.
    let eig_need = admit_dirichlet_eigs(24, 3, &ChebBudget::default()).expect("baseline admits");
    let mut eig_budget = ChebBudget::default();
    eig_budget.max_temp_bytes = eig_need.temp_bytes_admitted() - 1;
    assert!(matches!(
        admit_dirichlet_eigs(24, 3, &eig_budget),
        Err(ChebError::CapExceeded {
            what: "eigensolve temporary bytes",
            ..
        })
    ));
    eig_budget.max_temp_bytes = eig_need.temp_bytes_admitted();
    eig_budget.max_work_ops = eig_need.ops_admitted() - 1;
    assert!(matches!(
        admit_dirichlet_eigs(24, 3, &eig_budget),
        Err(ChebError::CapExceeded {
            what: "eigensolve work",
            ..
        })
    ));
    eig_budget.max_work_ops = eig_need.ops_admitted();
    admit_dirichlet_eigs(24, 3, &eig_budget).expect("exact eigensolve bounds admit");
    // Root scan: usize::MAX coefficients refuse before size arithmetic.
    assert!(matches!(
        admit_root_scan(usize::MAX, &ChebBudget::default()),
        Err(ChebError::CapExceeded { .. }) | Err(ChebError::Overflow { .. })
    ));
    let root_need = admit_root_scan(32, &ChebBudget::default()).expect("baseline admits");
    let mut root_budget = ChebBudget::default();
    root_budget.max_temp_bytes = root_need.temp_bytes_admitted() - 1;
    assert!(matches!(
        admit_root_scan(32, &root_budget),
        Err(ChebError::CapExceeded {
            what: "root-scan temporary bytes",
            ..
        })
    ));
    root_budget.max_temp_bytes = root_need.temp_bytes_admitted();
    root_budget.max_work_ops = root_need.ops_admitted() - 1;
    assert!(matches!(
        admit_root_scan(32, &root_budget),
        Err(ChebError::CapExceeded {
            what: "root-scan work",
            ..
        })
    ));
    root_budget.max_work_ops = root_need.ops_admitted();
    admit_root_scan(32, &root_budget).expect("exact root-scan bounds admit");
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

    let samples = core::cell::Cell::new(0usize);
    let refusal = with_cx(false, |cx| {
        try_build_budgeted(
            &|_| {
                samples.set(samples.get() + 1);
                1.0
            },
            0.0,
            1.0,
            17,
            Some(17),
            &ChebBudget::default(),
            cx,
        )
        .expect_err("rounded resume grids may not bypass max_degree")
    });
    assert!(matches!(refusal, ChebError::Shape { .. }), "{refusal}");
    assert_eq!(samples.get(), 0, "shape refusal precedes evaluation");
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

    let cancelled_during_sampling = with_gate_cx(|gate, cx| {
        try_build_budgeted(
            &|_| {
                gate.request();
                1.0
            },
            -1.0,
            1.0,
            64,
            None,
            &ChebBudget::default(),
            cx,
        )
        .expect("admitted sampling cancellation is a terminal state")
    });
    assert!(matches!(
        cancelled_during_sampling,
        BuildRun::Cancelled { .. }
    ));

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

    // Eigensolve: pre-cancelled drains with an EMPTY completed-estimate prefix.
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

        let exponent_range = Cheb1::from_coeffs(-1.0, 1.0, vec![f64::MAX, f64::from_bits(1)]);
        let refusal = exponent_range
            .roots_budgeted(&ChebBudget::default(), cx)
            .expect_err("lossy root normalization is a typed refusal, not an assertion panic");
        assert!(matches!(refusal, ChebError::Numerical { .. }), "{refusal}");
    });
}

/// cb-006 — local receipt determinism: identical budgeted runs produce
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

/// cb-007 — colleague admission + budgeted twin (slice 2): boundary
/// tables refuse before allocation; the budgeted path matches the
/// classic candidates bitwise; a pre-cancelled gate drains without
/// running the eigensolve.
#[test]
fn cb_007_colleague_budgeted() {
    use fs_cheb::colleague::ColleaguePolicy;
    // Shape/cap boundaries.
    assert!(matches!(
        fs_cheb::admit_colleague_roots(1, &ChebBudget::default()),
        Err(ChebError::Shape { .. })
    ));
    let mut tight = ChebBudget::default();
    tight.max_eigen_dim = 8;
    fs_cheb::admit_colleague_roots(9, &tight).expect("n = 8 admits at cap");
    assert!(matches!(
        fs_cheb::admit_colleague_roots(10, &tight),
        Err(ChebError::CapExceeded {
            what: "colleague matrix dimension",
            ..
        })
    ));
    assert!(matches!(
        fs_cheb::admit_colleague_roots(usize::MAX, &ChebBudget::default()),
        Err(ChebError::CapExceeded { .. }) | Err(ChebError::Overflow { .. })
    ));

    // Bitwise parity + cancellation drain.
    let poly = Cheb1::build(&|x: f64| (x - 0.2) * (x + 0.4) * (x - 0.9), -1.0, 1.0, 64);
    let classic = fs_cheb::colleague::colleague_roots(&poly, ColleaguePolicy::default());
    let budgeted = with_cx(false, |cx| {
        fs_cheb::colleague_roots_budgeted(
            &poly,
            ColleaguePolicy::default(),
            &ChebBudget::default(),
            cx,
        )
        .expect("admitted colleague run")
    });
    assert_eq!(budgeted.len(), classic.len());
    for (lhs, rhs) in budgeted.iter().zip(&classic) {
        assert_eq!(lhs.to_bits(), rhs.to_bits(), "candidate bit parity");
    }
    let refused = with_cx(true, |cx| {
        fs_cheb::colleague_roots_budgeted(
            &poly,
            ColleaguePolicy::default(),
            &ChebBudget::default(),
            cx,
        )
        .expect_err("pre-cancelled gate drains before the eigen tile")
    });
    assert!(matches!(refused, ChebError::Cancelled), "{refused}");
}

/// cb-008 — slice-2 admission tables for the remaining modules:
/// Cheb2 grids, Fourier synthesis, and the Orr–Sommerfeld eigensolve
/// all preflight their worst case with typed refusals (including the
/// panics the classic APIs would throw for shape violations).
#[test]
fn cb_008_remaining_module_admissions() {
    let b = ChebBudget::default();
    // Cheb2: domain/tol/rank shape refusals + grid caps.
    assert!(matches!(
        fs_cheb::admit_cheb2_build((0.0, 1.0, 0.0, 1.0), -1.0, 4, 64, &b),
        Err(ChebError::Shape { .. })
    ));
    assert!(matches!(
        fs_cheb::admit_cheb2_build((0.0, 1.0, 1.0, 1.0), 1e-12, 4, 64, &b),
        Err(ChebError::Domain { .. })
    ));
    assert!(matches!(
        fs_cheb::admit_cheb2_build((0.0, 1.0, 0.0, 1.0), 1e-12, 0, 64, &b),
        Err(ChebError::Shape { .. })
    ));
    let ok = fs_cheb::admit_cheb2_build((0.0, 1.0, 0.0, 1.0), 1e-12, 4, 64, &b)
        .expect("modest 2D build admits");
    assert_eq!(ok.samples_admitted(), 129 * 129, "ns = 128 grid side 129");
    assert!(matches!(
        fs_cheb::admit_cheb2_build((0.0, 1.0, 0.0, 1.0), 1e-12, 4, usize::MAX, &b),
        Err(ChebError::Overflow { .. }) | Err(ChebError::CapExceeded { .. })
    ));

    // Fourier: typed power-of-two refusal where the classic API panics.
    for bad in [0usize, 1, 3, 1000] {
        assert!(matches!(
            fs_cheb::admit_fourier_build(bad, &b),
            Err(ChebError::Shape { .. })
        ));
    }
    let four = fs_cheb::admit_fourier_build(1024, &b).expect("radix-2 admits");
    assert_eq!(four.samples_admitted(), 1024);
    let mut small = ChebBudget::default();
    small.max_samples = 512;
    assert!(matches!(
        fs_cheb::admit_fourier_build(1024, &small),
        Err(ChebError::CapExceeded {
            what: "Fourier samples",
            ..
        })
    ));

    // Orr–Sommerfeld: n >= 8, k in 1..=n, dimension/work caps.
    assert!(matches!(
        fs_cheb::admit_growth_rates(4, 1, &b),
        Err(ChebError::Shape { .. })
    ));
    assert!(matches!(
        fs_cheb::admit_growth_rates(48, 0, &b),
        Err(ChebError::Shape { .. })
    ));
    assert!(matches!(
        fs_cheb::admit_growth_rates(48, 49, &b),
        Err(ChebError::Shape { .. })
    ));
    fs_cheb::admit_growth_rates(48, 4, &b).expect("fixture-scale OS admits");
    assert!(matches!(
        fs_cheb::admit_growth_rates(usize::MAX, 1, &b),
        Err(ChebError::Overflow { .. }) | Err(ChebError::CapExceeded { .. })
    ));
}

/// cb-009 — slice-3 budgeted builder twins: the Fourier mirror is
/// bitwise-identical with typed refusals where the classic API panics;
/// the Cheb2 and Orr–Sommerfeld wrappers admit, drain at their tile
/// boundaries, and map solver failures to typed errors.
#[test]
fn cb_009_budgeted_builder_twins() {
    with_cx(false, |cx| {
        // Fourier: bitwise mirror parity.
        let wave = |theta: f64| fs_math::det::sin(theta) + 0.25 * fs_math::det::cos(3.0 * theta);
        let classic = fs_cheb::fourier::FourierSeries::build(&wave, 64);
        let budgeted = fs_cheb::fourier_build_budgeted(&wave, 64, &ChebBudget::default(), cx)
            .expect("admitted synthesis");
        assert_eq!(budgeted.n, classic.n);
        assert_eq!(budgeted.coeffs.len(), classic.coeffs.len());
        for (lhs, rhs) in budgeted.coeffs.iter().zip(&classic.coeffs) {
            assert_eq!(lhs.re.to_bits(), rhs.re.to_bits(), "spectrum bit parity");
            assert_eq!(lhs.im.to_bits(), rhs.im.to_bits(), "spectrum bit parity");
        }
        // Typed refusals where the classic API panics.
        assert!(matches!(
            fs_cheb::fourier_build_budgeted(&wave, 100, &ChebBudget::default(), cx),
            Err(ChebError::Shape { .. })
        ));
        // theta = 0 IS the k = 0 sample point, so the pole is hit
        // exactly and the refusal is deterministic.
        let singular = |theta: f64| 1.0 / theta;
        assert!(matches!(
            fs_cheb::fourier_build_budgeted(&singular, 64, &ChebBudget::default(), cx),
            Err(ChebError::NonFinite { .. })
        ));

        // Cheb2: admitted wrapper matches the classic factors bitwise.
        let f2 = |x: f64, y: f64| fs_math::det::sin(2.0 * x) * fs_math::det::cos(y);
        let classic2 = fs_cheb::cheb2::Cheb2::build(&f2, (0.0, 1.0, 0.0, 2.0), 1e-10, 6, 64);
        let budgeted2 = fs_cheb::cheb2_build_budgeted(
            &f2,
            (0.0, 1.0, 0.0, 2.0),
            1e-10,
            6,
            64,
            &ChebBudget::default(),
            cx,
        )
        .expect("admitted 2D build");
        assert_eq!(budgeted2.rank(), classic2.rank());
        assert_eq!(
            budgeted2.residual.to_bits(),
            classic2.residual.to_bits(),
            "accuracy ledger bit parity"
        );

        // Orr–Sommerfeld: parity with the classic solver + typed
        // parameter refusal.
        let classic_os =
            fs_cheb::orr_sommerfeld::growth_rates(5772.0, 1.02, 32, 3).expect("classic OS");
        let budgeted_os =
            fs_cheb::growth_rates_budgeted(5772.0, 1.02, 32, 3, &ChebBudget::default(), cx)
                .expect("admitted OS");
        assert_eq!(budgeted_os.len(), classic_os.len());
        for (lhs, rhs) in budgeted_os.iter().zip(&classic_os) {
            assert_eq!(lhs.re.to_bits(), rhs.re.to_bits(), "growth-rate bit parity");
            assert_eq!(lhs.im.to_bits(), rhs.im.to_bits(), "growth-rate bit parity");
        }
        assert!(matches!(
            fs_cheb::growth_rates_budgeted(f64::NAN, 1.02, 32, 3, &ChebBudget::default(), cx),
            Err(ChebError::NonFinite { .. })
        ));
    });

    // Cancellation drains for all three twins (pre-cancelled gates).
    let wave = |theta: f64| fs_math::det::sin(theta);
    let f2 = |x: f64, y: f64| x + y;
    assert!(matches!(
        with_cx(true, |cx| fs_cheb::fourier_build_budgeted(
            &wave,
            64,
            &ChebBudget::default(),
            cx
        )),
        Err(ChebError::Cancelled)
    ));
    assert!(matches!(
        with_cx(true, |cx| fs_cheb::cheb2_build_budgeted(
            &f2,
            (0.0, 1.0, 0.0, 1.0),
            1e-10,
            4,
            64,
            &ChebBudget::default(),
            cx
        )),
        Err(ChebError::Cancelled)
    ));
    assert!(matches!(
        with_cx(true, |cx| fs_cheb::growth_rates_budgeted(
            5772.0,
            1.02,
            32,
            3,
            &ChebBudget::default(),
            cx
        )),
        Err(ChebError::Cancelled)
    ));
}

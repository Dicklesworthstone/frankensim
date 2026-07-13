//! Abstraction-ladder conformance (the knh1.4 bead; runs under
//! `abstraction-ladder`; colors re-graded by bead y6yv). Acceptance: RB
//! estimators empirically contain the observed errors on the canonical
//! parametric elliptic fixture (with an explicit roundoff slack) — and
//! travel as ESTIMATED dispersion, never Verified, because they are
//! f64-evaluated without outward rounding or an independent oracle;
//! the leak alarm never lets an operator act on an abstraction whose
//! bound exceeds tolerance — it auto-drills; every invalid input is a
//! structured refusal, never a panic or a colored NaN; the kill
//! measurement (RB coverage) is ledgered above the 20% beachhead
//! floor; queries replay bit-equal (G5).
#![cfg(feature = "abstraction-ladder")]

use core::mem::size_of;
use fs_alloc::OperationMemoryLease;
use fs_evidence::{Color, NumericalKind, validate_color_payload};
use fs_exec::{Budget, CancelGate, RunId, TilePool};
use fs_surrogate::ladder::{
    ConceptLevel, Ladder, MAX_CONCEPT_POINTS, MAX_COVERAGE_AXIS_POINTS, MAX_COVERAGE_QUERIES,
    MAX_RB_DIM, MAX_RB_LEVELS, MAX_SYNCHRONOUS_COVERAGE_WORK_UNITS, MAX_TRAINING_BYTES,
    MAX_TRAINING_WORK_UNITS, RbCoverageOutcome, RbCoveragePlan, RbLevel, SurrogateError,
    TruthModel, rb_coverage, rb_coverage_scoped,
};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-surrogate/ladder\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

#[test]
fn la_001_rb_estimators_empirically_contain_fixture_errors() {
    // Empirical fixture check: |s_truth − s_rb| ≤ qoi_estimator plus an
    // explicit absolute roundoff slack across a μ battery. Both sides use
    // this crate's f64 path, so this is useful regression evidence, not an
    // independent enclosure or a certifier test (y6yv).
    const ABSOLUTE_ROUNDOFF_SLACK: f64 = 1e-14;
    let truth = TruthModel::new(200).expect("bounded dimension");
    for k in [2usize, 4, 6] {
        let rb = RbLevel::train(&truth, (0.0, 4.0), k).expect("canonical training");
        let mut worst_ratio = 0.0f64;
        for i in 0..25 {
            let mu = 4.0 * f64::from(i) / 24.0;
            let (u_rb, s_rb, _, qoi_bound) = rb.query(mu).expect("in-range query");
            let galerkin_defect = (s_rb
                - truth
                    .energy(&u_rb, &u_rb, mu)
                    .expect("finite reduced energy"))
            .abs();
            let s_true = truth
                .compliance(&truth.solve(mu).expect("coercive solve"))
                .expect("finite compliance");
            let err = (s_true - s_rb).abs();
            assert!(
                err <= qoi_bound + ABSOLUTE_ROUNDOFF_SLACK,
                "k={k}, mu={mu}: |{s_true} - {s_rb}| = {err:.3e} exceeds estimator \
                 {qoi_bound:.3e} plus slack {ABSOLUTE_ROUNDOFF_SLACK:.1e}"
            );
            assert!(
                qoi_bound >= galerkin_defect,
                "the compliance estimator must retain floating Galerkin defect"
            );
            if qoi_bound > 1e-15 {
                worst_ratio = worst_ratio.max(err / qoi_bound);
            }
        }
        println!(
            "{{\"metric\":\"rb-bound\",\"k\":{k},\"worst_effectivity_inverse\":{worst_ratio:.4}}}"
        );
    }
    // For this one fixed off-training probe the independently trained fixture
    // spaces happen to tighten at each adjacent basis size. These spaces are
    // not nested, so this observation is deliberately not a general monotonicity
    // claim.
    let bounds: Vec<f64> = [2usize, 4, 6]
        .iter()
        .map(|&k| {
            RbLevel::train(&truth, (0.0, 4.0), k)
                .expect("training")
                .query(1.7)
                .expect("query")
                .3
        })
        .collect();
    // The estimator improved to the point where k=4 and k=6 both land at
    // the roundoff floor (~1e-16); strict ordering BETWEEN two at-floor
    // values is noise, not tightening. The observation keeps its teeth:
    // every adjacent pair must tighten OR both sit at the floor, and the
    // first refinement must still tighten visibly from above the floor.
    const ROUNDOFF_FLOOR: f64 = 1e-15; // matches the effectivity guard above
    assert!(
        bounds
            .windows(2)
            .all(|pair| pair[1] < pair[0]
                || (pair[0] <= ROUNDOFF_FLOOR && pair[1] <= ROUNDOFF_FLOOR)),
        "the canonical off-grid fixture stopped showing its observed adjacent tightening \
         (or floor containment): {bounds:?}"
    );
    assert!(
        bounds[0] > ROUNDOFF_FLOOR && bounds[1] < bounds[0],
        "the first refinement must still tighten from above the floor: {bounds:?}"
    );
    verdict(
        "la-001",
        "the f64 RB QoI estimator empirically contains the canonical fixture errors with \
         1e-14 absolute slack for k in {2,4,6}; one off-grid probe shows adjacent \
         tightening, without claiming nested-space monotonicity or certification",
    );
}

#[test]
fn la_002_leak_alarm_never_lets_a_leak_answer() {
    // THE PROPERTY: whatever rung you start at, the answer's quantified
    // width meets the tolerance or came from level 0 — a leaking rung
    // NEVER answers; it is recorded and descended past.
    let ladder = Ladder::build(200, (0.0, 4.0), &[6, 2], true).expect("canonical ladder");
    for i in 0..15 {
        let mu = 0.2 + 3.6 * f64::from(i) / 14.0;
        for tol in [1e-3, 1e-6, 1e-10, 1e-14] {
            let ans = ladder
                .at_level(ladder.top())
                .expect("declared top level")
                .query(mu, tol)
                .expect("valid query");
            let Color::Estimated { dispersion, .. } = ans.color() else {
                assert!(
                    matches!(ans.color(), Color::Estimated { .. }),
                    "y6yv: the ladder only mints Estimated, got {:?}",
                    ans.color()
                );
                continue;
            };
            assert!(
                *dispersion <= tol || ans.level_used() == 0,
                "mu={mu}, tol={tol:.0e}: dispersion {dispersion:.2e} from level {}",
                ans.level_used()
            );
            for &l in ans.leaks() {
                assert!(
                    l > ans.level_used(),
                    "leak {l} above answer {}",
                    ans.level_used()
                );
            }
        }
    }
    // A tolerance below every rung's achievable bound forces full
    // descent with the whole ordered leak trail.
    let ans = ladder
        .at_level(ladder.top())
        .expect("declared top level")
        .query(1.3, 1e-32)
        .expect("valid query");
    assert_eq!(ans.level_used(), 0, "ultra-tight tol reaches the truth");
    assert_eq!(ans.leaks(), [3, 2, 1], "every rung leaked, in order");
    verdict(
        "la-002",
        "across 15 mu x 4 tolerances no leaking rung ever answers: the quantified \
         dispersion meets tol or level 0 answered; the leak trail is complete and ordered",
    );
}

#[test]
fn la_003_no_color_masquerades_as_verified() {
    let ladder = Ladder::build(200, (0.0, 4.0), &[6, 2], true).expect("canonical ladder");
    // A loose query the concept rung CAN answer: honest Estimated.
    let loose = ladder
        .at_level(ladder.top())
        .expect("declared top level")
        .query(2.0, 0.5)
        .expect("valid query");
    assert_eq!(
        loose.level_used(),
        ladder.top(),
        "the concept rung answered"
    );
    assert!(
        matches!(loose.color(), Color::Estimated { estimator, .. }
            if estimator == "fs-surrogate.concept-cross-rung-v1"),
        "the concept answer is estimated, calibrated by cross-rung probes: {:?}",
        loose.color()
    );
    // The same query demanding tight accuracy skips the concept rung
    // and answers from an RB rung — whose color is ALSO Estimated
    // (y6yv: an f64-evaluated textbook bound is not an enclosure), with
    // the estimator naming exactly what it is.
    let strict = ladder
        .at_level(ladder.top())
        .expect("declared top level")
        .query(2.0, 1e-6)
        .expect("valid query");
    assert!(
        matches!(strict.color(), Color::Estimated { estimator, .. }
            if estimator == "fs-surrogate.rb-a-posteriori-f64-v1"),
        "the RB answer names its authority honestly: {:?}",
        strict.color()
    );
    assert!(
        strict.leaks().contains(&ladder.top()),
        "the concept rung is recorded as leaked, not silently skipped"
    );
    // Level 0 names DECLARED truth semantics and makes no spread claim:
    // an unproved floating solve must carry unbounded dispersion.
    let full = ladder
        .at_level(0)
        .expect("truth level")
        .query(2.0, 1e-32)
        .expect("valid query");
    assert!(
        matches!(full.color(), Color::Estimated { estimator, dispersion }
            if estimator == "fs-surrogate.fe-truth-f64-no-certificate-v1"
                && dispersion.is_infinite()),
        "level 0 is a declared-semantics point value, not a verified interval: {:?}",
        full.color()
    );
    for answer in [&loose, &strict, &full] {
        validate_color_payload(answer.color()).expect("ladder emits valid color payloads");
    }
    verdict(
        "la-003",
        "every rung's color is Estimated with an estimator naming its true authority — \
         nothing in the ladder can mint Verified until outward-rounded certificates exist",
    );
}

#[test]
fn la_004_kill_measurement_rb_coverage() {
    let ladder = Ladder::build(200, (0.0, 4.0), &[6, 2], false).expect("canonical ladder");
    let mus: Vec<f64> = (0..12).map(|i| 4.0 * f64::from(i) / 11.0).collect();
    let tols = [1e-2, 1e-4, 1e-6, 1e-8];
    let coverage = rb_coverage(&ladder, &mus, &tols).expect("valid battery");
    let reference_view = ladder
        .at_level(ladder.rb_level_count())
        .expect("coarsest RB level");
    let mut reference_covered = 0usize;
    for &mu in &mus {
        for &tol in &tols {
            let answer = reference_view
                .query(mu, tol)
                .expect("reference LevelView battery query");
            if answer.level_used() >= 1 {
                reference_covered += 1;
            }
        }
    }
    #[allow(clippy::cast_precision_loss)]
    let reference_coverage = reference_covered as f64 / (mus.len() * tols.len()) as f64;
    assert_eq!(
        coverage.to_bits(),
        reference_coverage.to_bits(),
        "one-solve-per-mu classification preserves LevelView coverage semantics"
    );
    println!(
        "{{\"metric\":\"kill-measurement\",\"rb_coverage\":{coverage:.3},\
         \"floor\":0.2,\"battery\":\"12 mu x 4 tol\"}}"
    );
    assert!(
        coverage >= 0.2,
        "the beachhead clears the kill floor: {coverage}"
    );
    verdict(
        "la-004",
        "RB coverage of the 48-point query battery ledgered and above the 0.2 kill floor",
    );
}

#[test]
fn la_004b_scoped_coverage_replays_across_worker_counts() {
    let ladder = Ladder::build(200, (0.0, 4.0), &[6, 2], false).expect("canonical ladder");
    let mus: Vec<f64> = (0..12).map(|i| 4.0 * f64::from(i) / 11.0).collect();
    let tolerances = [1e-2, 1e-4, 1e-6, 1e-8];
    let synchronous = rb_coverage(&ladder, &mus, &tolerances).expect("small oracle battery");
    let plan = RbCoveragePlan::new(&ladder, &mus, &tolerances).expect("exact typed plan");

    let mut observed = Vec::new();
    for workers in [1, 4] {
        let pool = TilePool::for_host(workers, 0xC0_4E_2A_6E);
        let cx = asupersync::Cx::for_testing();
        let gate = CancelGate::new();
        let lease = OperationMemoryLease::bounded(64 * 1024 * 1024);
        let run = rb_coverage_scoped(
            &cx,
            &pool,
            &gate,
            RunId(0x4B),
            Budget::INFINITE,
            &lease,
            &ladder,
            &plan,
        )
        .expect("scoped coverage completes");
        let RbCoverageOutcome::Complete {
            coverage,
            covered_queries,
            authority,
            receipt,
        } = run.outcome()
        else {
            panic!("healthy scoped coverage must finalize: {:?}", run.outcome());
        };
        assert_eq!(coverage.to_bits(), synchronous.to_bits());
        #[allow(clippy::cast_precision_loss)]
        let counted_coverage = *covered_queries as f64 / plan.total_queries() as f64;
        assert_eq!(
            counted_coverage.to_bits(),
            synchronous.to_bits(),
            "the complete outcome count and fraction agree exactly"
        );
        assert!(matches!(
            authority,
            Color::Estimated {
                estimator,
                dispersion
            } if estimator == "fs-surrogate.rb-coverage-f64-v2" && dispersion.is_infinite()
        ));
        assert!(receipt.is_finalized());
        assert_eq!(receipt.completed_parameter_prefix(), mus.len());
        assert!(receipt.current_parameter().is_none());
        assert_eq!(receipt.plan(), &plan);
        assert_eq!(run.report().completed, mus.len() as u64);
        assert_eq!(run.report().total, mus.len() as u64);
        assert_eq!(run.memory_receipt().used_bytes, 0);
        assert_eq!(run.memory_receipt().refusals, 0);
        assert!(pool.arena_pool().stats().quiescent());
        observed.push((coverage.to_bits(), *covered_queries));
    }
    assert_eq!(observed[0], observed[1]);
    verdict(
        "la-004b",
        "production coverage is bit-identical to the synchronous oracle across worker counts; \
         only measured scheduling diagnostics vary",
    );
}

#[test]
fn la_004c_g4_cancellation_drains_retains_prefix_and_reuses_pool() {
    let ladder = Ladder::build(2_048, (0.0, 4.0), &[4, 2], false).expect("canonical ladder");
    let mus: Vec<f64> = (0..8).map(|i| 4.0 * f64::from(i) / 7.0).collect();
    let tolerances = [1e-4, 1e-6, 1e-8];
    let plan = RbCoveragePlan::new(&ladder, &mus, &tolerances).expect("exact typed plan");
    let pool = TilePool::for_host(1, 0xC4_4C_E1);

    let pre_cancelled = asupersync::Cx::for_testing();
    pre_cancelled.set_cancel_requested(true);
    let pre_gate = CancelGate::new();
    let pre_lease = OperationMemoryLease::bounded(64 * 1024 * 1024);
    let pre = rb_coverage_scoped(
        &pre_cancelled,
        &pool,
        &pre_gate,
        RunId(0x4C),
        Budget::INFINITE,
        &pre_lease,
        &ladder,
        &plan,
    )
    .expect("pre-cancellation is an incomplete outcome");
    let RbCoverageOutcome::Incomplete { no_claim, receipt } = pre.outcome() else {
        panic!("pre-cancelled work cannot publish coverage");
    };
    assert_eq!(no_claim.kind, NumericalKind::NoClaim);
    assert_eq!(receipt.completed_parameter_prefix(), 0);
    assert!(!receipt.is_finalized());
    assert_eq!(pre.memory_receipt().used_bytes, 0);
    assert!(pool.arena_pool().stats().quiescent());

    let mut replay_receipt = None;
    let mut observed_latency = None;
    for poll_quota in [8, 16, 32, 64] {
        let cx = asupersync::Cx::for_testing_with_budget(Budget {
            poll_quota,
            ..Budget::INFINITE
        });
        let gate = CancelGate::new();
        let lease = OperationMemoryLease::bounded(64 * 1024 * 1024);
        let run = rb_coverage_scoped(
            &cx,
            &pool,
            &gate,
            RunId(0x4C),
            Budget::INFINITE,
            &lease,
            &ladder,
            &plan,
        )
        .expect("budget exhaustion is a drained incomplete outcome");
        let RbCoverageOutcome::Incomplete { no_claim, receipt } = run.outcome() else {
            panic!("bounded cancellation storm unexpectedly published coverage");
        };
        assert_eq!(no_claim.kind, NumericalKind::NoClaim);
        assert!(!receipt.is_finalized());
        assert!(receipt.completed_parameter_prefix() < mus.len());
        assert_eq!(run.memory_receipt().used_bytes, 0);
        assert!(pool.arena_pool().stats().quiescent());
        observed_latency = observed_latency.or_else(|| run.report().cancel_latency_p99_ns());
        if poll_quota == 32 {
            replay_receipt = Some(receipt.clone());
        }
    }
    assert!(
        observed_latency.is_some_and(|latency| latency < 5_000_000_000),
        "a mid-kernel cancellation must produce a finite drained latency sample"
    );

    let replay_cx = asupersync::Cx::for_testing_with_budget(Budget {
        poll_quota: 32,
        ..Budget::INFINITE
    });
    let replay_gate = CancelGate::new();
    let replay_lease = OperationMemoryLease::bounded(64 * 1024 * 1024);
    let replay = rb_coverage_scoped(
        &replay_cx,
        &pool,
        &replay_gate,
        RunId(0x4C),
        Budget::INFINITE,
        &replay_lease,
        &ladder,
        &plan,
    )
    .expect("deterministic cancellation replay");
    assert_eq!(
        replay.outcome().receipt(),
        replay_receipt.as_ref().expect("quota-32 receipt retained"),
        "same logical run and poll budget retain the same canonical prefix/progress"
    );

    let healthy_cx = asupersync::Cx::for_testing();
    let healthy_gate = CancelGate::new();
    let healthy_lease = OperationMemoryLease::bounded(64 * 1024 * 1024);
    let healthy = rb_coverage_scoped(
        &healthy_cx,
        &pool,
        &healthy_gate,
        RunId(0x4C),
        Budget::INFINITE,
        &healthy_lease,
        &ladder,
        &plan,
    )
    .expect("pool remains reusable after cancellation storms");
    assert!(matches!(
        healthy.outcome(),
        RbCoverageOutcome::Complete { .. }
    ));
    assert!(pool.arena_pool().stats().quiescent());
    verdict(
        "la-004c",
        "G4 pre-cancel and budget storms drain with no fraction, retain a replayable prefix, \
         release memory, and leave the pool reusable",
    );
}

#[test]
fn la_005_g5_determinism_and_invisibility() {
    let ladder = Ladder::build(150, (0.0, 4.0), &[5], true).expect("canonical ladder");
    let level = ladder.at_level(1).expect("declared RB level");
    let a = level.query(2.3, 1e-6).expect("valid query");
    let b = level.query(2.3, 1e-6).expect("valid query");
    assert!(a.value().to_bits() == b.value().to_bits() && a.level_used() == b.level_used());
    let quiet = level.query(2.3, 1e-2).expect("valid query");
    assert_eq!(quiet.level_used(), 1);
    assert!(quiet.leaks().is_empty(), "no leak, no descent, no noise");
    verdict(
        "la-005",
        "queries replay bit-equal (G5); a satisfiable query answers at its rung with an \
         empty leak trail — the ladder is invisible until it leaks",
    );
}

#[test]
#[allow(clippy::too_many_lines)] // one hostile public-API refusal matrix
fn la_006_adversarial_inputs_are_structured_refusals() {
    assert!(matches!(
        TruthModel::new(0),
        Err(SurrogateError::InvalidDimension { n: 0 })
    ));
    assert!(matches!(
        TruthModel::new(usize::MAX),
        Err(SurrogateError::InvalidDimension { .. })
    ));
    let truth = TruthModel::new(50).expect("bounded");
    for bad in [-1.0, -2.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            truth.solve(bad),
            Err(SurrogateError::InvalidCoercivity { .. })
        ));
    }
    for bad_range in [(4.0, 0.0), (1.0, 1.0), (f64::NAN, 2.0)] {
        assert!(matches!(
            RbLevel::train(&truth, bad_range, 3),
            Err(SurrogateError::InvalidRange { .. })
        ));
    }
    assert!(matches!(
        RbLevel::train(&truth, (-3.0, 2.0), 3),
        Err(SurrogateError::InvalidCoercivity { .. })
    ));

    let finite = vec![1.0; truth.n()];
    assert!(matches!(
        truth.energy(&[], &finite, 0.0),
        Err(SurrogateError::InvalidVectorLength { .. })
    ));
    assert!(matches!(
        truth.compliance(&[]),
        Err(SurrogateError::InvalidVectorLength { .. })
    ));
    let mut non_finite = finite.clone();
    non_finite[7] = f64::NAN;
    assert!(matches!(
        truth.energy(&non_finite, &finite, 0.0),
        Err(SurrogateError::NonFiniteInput { index: 7, .. })
    ));
    assert!(matches!(
        truth.compliance(&non_finite),
        Err(SurrogateError::NonFiniteInput { index: 7, .. })
    ));
    let overflowing = vec![f64::MAX; truth.n()];
    assert!(matches!(
        truth.energy(&overflowing, &overflowing, 0.0),
        Err(SurrogateError::NonFiniteDerived { .. })
    ));
    assert!(matches!(
        truth.compliance(&overflowing),
        Err(SurrogateError::NonFiniteDerived { .. })
    ));
    assert!(matches!(
        truth.solve(f64::MAX),
        Err(SurrogateError::NonFiniteDerived { .. })
    ));

    let rb = RbLevel::train(&truth, (0.0, 4.0), 3).expect("training");
    assert!(matches!(
        rb.query(5.0),
        Err(SurrogateError::OutOfRange { .. })
    ));
    assert!(matches!(
        rb.query(f64::NAN),
        Err(SurrogateError::InvalidCoercivity { .. })
    ));
    let concept = ConceptLevel::train(&rb, 5, 7).expect("concept training");
    assert!(matches!(
        concept.lookup(-0.25),
        Err(SurrogateError::OutOfRange { .. })
    ));
    assert!(matches!(
        concept.lookup(f64::NAN),
        Err(SurrogateError::InvalidCoercivity { .. })
    ));
    assert!(matches!(
        ConceptLevel::train(&rb, 1, 7),
        Err(SurrogateError::InvalidGeometry { .. })
    ));
    assert!(matches!(
        ConceptLevel::train(&rb, 5, 0),
        Err(SurrogateError::InvalidGeometry { .. })
    ));

    let ladder = Ladder::build(50, (0.0, 4.0), &[3], false).expect("ladder");
    assert!(matches!(
        ladder.at_level(ladder.top() + 1),
        Err(SurrogateError::InvalidLevel { .. })
    ));
    let level = ladder.at_level(1).expect("declared level");
    for bad_tol in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            level.query(2.0, bad_tol),
            Err(SurrogateError::InvalidTolerance { .. })
        ));
    }
    assert!(matches!(
        level.query(9.0, 1e-4),
        Err(SurrogateError::OutOfRange { .. })
    ));
    verdict(
        "la-006",
        "public truth, RB, concept, and level APIs refuse bad shapes, non-finite arithmetic, \
         non-coercive and out-of-range parameters, invalid tolerances, and invalid levels",
    );
}

#[test]
#[allow(clippy::too_many_lines)] // one preflight matrix spanning all coupled resource invariants
fn la_007_training_preflight_is_bounded_and_representable() {
    let truth = TruthModel::new(1_024).expect("bounded truth");
    for dimension in [0, MAX_RB_DIM + 1] {
        assert!(matches!(
            RbLevel::train(&truth, (0.0, 4.0), dimension),
            Err(SurrogateError::InvalidRbDimension { .. })
        ));
    }
    let small_truth = TruthModel::new(4).expect("small truth");
    assert!(matches!(
        RbLevel::train(&small_truth, (0.0, 4.0), 5),
        Err(SurrogateError::InvalidRbDimension { maximum: 4, .. })
    ));

    assert!(matches!(
        Ladder::build(64, (0.0, 4.0), &[], false),
        Err(SurrogateError::InvalidGeometry { .. })
    ));
    let too_many = (1..=MAX_RB_LEVELS + 1).rev().collect::<Vec<_>>();
    assert!(matches!(
        Ladder::build(64, (0.0, 4.0), &too_many, false),
        Err(SurrogateError::TooManyRbLevels { .. })
    ));
    for dimensions in [&[4, 4][..], &[2, 4][..]] {
        assert!(matches!(
            Ladder::build(64, (0.0, 4.0), dimensions, false),
            Err(SurrogateError::NonDecreasingFidelity { .. })
        ));
    }

    let memory_n = MAX_TRAINING_BYTES / (MAX_RB_DIM * size_of::<f64>()) + 1;
    assert!(matches!(
        Ladder::build(memory_n, (0.0, 4.0), &[MAX_RB_DIM], false),
        Err(SurrogateError::BudgetExceeded {
            resource: "aggregate ladder training memory",
            ..
        })
    ));
    let work_dimension = 64;
    let work_n = MAX_TRAINING_WORK_UNITS / (work_dimension * work_dimension) + 1;
    assert!(
        work_n * (work_dimension + 6) * size_of::<f64>()
            + work_dimension * work_dimension * size_of::<f64>()
            <= MAX_TRAINING_BYTES,
        "work hostile fixture must clear the independent memory cap"
    );
    assert!(matches!(
        Ladder::build(work_n, (0.0, 4.0), &[work_dimension], false),
        Err(SurrogateError::BudgetExceeded {
            resource: "aggregate ladder training work",
            ..
        })
    ));

    let narrow = (1.0, f64::from_bits(1.0f64.to_bits() + 8));
    assert!(matches!(
        RbLevel::train(&truth, narrow, 16),
        Err(SurrogateError::UnrepresentableGrid { .. })
    ));
    let narrow_rb = RbLevel::train(&truth, narrow, 2).expect("two distinct endpoints");
    assert_eq!(narrow_rb.dim(), 1, "narrow snapshots collapse to one mode");
    assert!(matches!(
        Ladder::build(truth.n(), narrow, &[2, 1], false),
        Err(SurrogateError::NonDecreasingRetainedFidelity {
            previous: 1,
            current: 1,
            ..
        })
    ));
    assert!(matches!(
        ConceptLevel::train(&narrow_rb, MAX_CONCEPT_POINTS, 1),
        Err(SurrogateError::UnrepresentableGrid { .. })
    ));
    let large_n = MAX_TRAINING_WORK_UNITS / (2 * MAX_CONCEPT_POINTS) + 1;
    let large_truth = TruthModel::new(large_n).expect("bounded large truth");
    let large_rb = RbLevel::train(&large_truth, (0.0, 4.0), 1).expect("bounded large RB training");
    assert!(matches!(
        ConceptLevel::train(&large_rb, MAX_CONCEPT_POINTS, MAX_CONCEPT_POINTS),
        Err(SurrogateError::BudgetExceeded {
            resource: "concept training work",
            ..
        })
    ));
    let adjacent = (1.0, f64::from_bits(1.0f64.to_bits() + 1));
    assert!(matches!(
        Ladder::build(64, adjacent, &[2], false),
        Err(SurrogateError::UnrepresentableGrid { .. })
    ));

    let ladder = Ladder::build(64, (0.0, 4.0), &[4, 2], true).expect("valid ladder");
    assert!(
        ladder
            .rb_levels()
            .windows(2)
            .all(|pair| pair[0].dim() > pair[1].dim()),
        "stored retained dimensions are strictly decreasing"
    );
    for level in ladder.rb_levels() {
        assert_eq!(level.family(), ladder.family());
    }
    assert_eq!(
        ladder.concept().expect("concept rung").family(),
        ladder.family()
    );
    verdict(
        "la-007",
        "rung order/count/dimensions and aggregate memory/work are preflighted before training; \
         retained dimensions must also descend; unrepresentable grids refuse; every sealed rung \
         shares one family/range identity",
    );
}

#[test]
fn la_008_coverage_batteries_are_bounded_and_nonempty() {
    let ladder = Ladder::build(64, (0.0, 4.0), &[4], false).expect("valid ladder");
    assert!(matches!(
        rb_coverage(&ladder, &[], &[1e-4]),
        Err(SurrogateError::EmptyCoverageBattery { axis: "mu" })
    ));
    assert!(matches!(
        rb_coverage(&ladder, &[0.0], &[]),
        Err(SurrogateError::EmptyCoverageBattery { axis: "tolerance" })
    ));
    let oversized_axis = vec![0.0; MAX_COVERAGE_AXIS_POINTS + 1];
    assert!(matches!(
        rb_coverage(&ladder, &oversized_axis, &[1e-4]),
        Err(SurrogateError::CoverageAxisTooLarge { axis: "mu", .. })
    ));
    let mus = vec![0.0; 1_001];
    let tols = vec![1e-4; 1_000];
    assert!(mus.len() * tols.len() > MAX_COVERAGE_QUERIES);
    assert!(matches!(
        rb_coverage(&ladder, &mus, &tols),
        Err(SurrogateError::CoverageProductTooLarge { .. })
    ));
    assert!(matches!(
        rb_coverage(&ladder, &[9.0], &[1e-4]),
        Err(SurrogateError::OutOfRange { .. })
    ));
    for bad_tol in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            rb_coverage(&ladder, &[2.0], &[bad_tol]),
            Err(SurrogateError::InvalidTolerance { .. })
        ));
    }
    verdict(
        "la-008",
        "coverage refuses empty, oversized, over-product, out-of-range, and invalid-tolerance \
         batteries before executing the Cartesian query workload",
    );
}

#[test]
fn la_009_concept_dispersion_inherits_lower_rung_uncertainty() {
    let truth = TruthModel::new(200).expect("bounded truth");
    let rb = RbLevel::train(&truth, (0.0, 4.0), 1).expect("one-mode RB");
    let concept = ConceptLevel::train(&rb, 2, 1).expect("one-probe concept");
    let mu = 2.0;
    let concept_value = concept.lookup(mu).expect("in-range concept query");
    let (_, rb_value, _, rb_qoi_estimator) = rb.query(mu).expect("in-range RB query");
    let cross_rung_only = (concept_value - rb_value).abs();
    let expected_total = cross_rung_only + rb_qoi_estimator;

    assert!(
        rb_qoi_estimator > 0.0,
        "the hostile one-mode lower rung must retain visible uncertainty"
    );
    assert_eq!(
        concept.dispersion().to_bits(),
        expected_total.to_bits(),
        "one probe stores |concept-RB| plus the lower RB QoI estimator"
    );
    assert!(
        concept.dispersion() > cross_rung_only,
        "matching an inaccurate lower rung cannot erase its uncertainty"
    );

    let ladder = Ladder::build(200, (0.0, 4.0), &[1], true).expect("concept ladder");
    let query_mu = 1.3;
    let concept = ladder.concept().expect("concept rung");
    let concept_value = concept.lookup(query_mu).expect("off-probe lookup");
    let (_, rb_value, _, rb_qoi_estimator) = ladder.rb_levels()[0]
        .query(query_mu)
        .expect("query-local lower rung");
    let expected_query_dispersion = concept
        .dispersion()
        .max((concept_value - rb_value).abs() + rb_qoi_estimator);
    let answer = ladder
        .at_level(ladder.top())
        .expect("concept level")
        .query(query_mu, f64::MAX)
        .expect("finite off-probe query");
    assert_eq!(answer.level_used(), ladder.top());
    assert!(matches!(
        answer.color(),
        Color::Estimated { dispersion, .. }
            if dispersion.to_bits() == expected_query_dispersion.to_bits()
    ));
    verdict(
        "la-009",
        "concept dispersion includes calibrated and query-local cross-rung discrepancy plus \
         lower-rung QoI uncertainty, so an inaccurate RB cannot bypass descent",
    );
}

#[test]
fn la_010_coverage_reuses_per_mu_solves_and_caps_aggregate_work() {
    let ladder = Ladder::build(16_384, (0.0, 4.0), &[2, 1], false).expect("bounded ladder");

    // One parameter with the maximum tolerance axis is admitted because RB
    // solves are performed once for that parameter, not once per tolerance.
    let tolerances = vec![1.0e-4; MAX_COVERAGE_AXIS_POINTS];
    let reused = rb_coverage(&ladder, &[2.0], &tolerances)
        .expect("tolerance classification reuses one per-mu descent");
    assert!(reused.is_finite());

    // The orthogonal hostile shape remains bounded in the solver dimension:
    // many parameters require distinct solves even with one tolerance.
    let mus = vec![2.0; MAX_COVERAGE_AXIS_POINTS];
    assert!(matches!(
        rb_coverage(&ladder, &mus, &[1.0e-4]),
        Err(SurrogateError::BudgetExceeded {
            resource: "synchronous RB coverage total work",
            limit: MAX_SYNCHRONOUS_COVERAGE_WORK_UNITS,
            ..
        })
    ));
    verdict(
        "la-010",
        "coverage reuses one bounded descent across every tolerance for each parameter, while \
         refusing excessive aggregate work across independently solved parameters",
    );
}

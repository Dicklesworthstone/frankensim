//! The engine's own battery (bead frankensim-4nh8): determinism of the
//! case stream, shrink-lattice sanity, the seeded-violation self-test
//! (a planted law break must shrink to its known minimal kernel), and
//! replay-path equivalence. JSONL verdicts per house style.

use fs_propcheck::{
    MinimizeBudget, MinimizeError, Shrink, Stream, check, minimize, minimize_with_budget,
};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-propcheck\",\"case\":\"{case}\",\"verdict\":\"pass\",\"detail\":\"{detail}\"}}"
    );
}

#[test]
fn case_streams_are_deterministic_and_decorrelated() {
    // Same (seed, index) twice: identical draws.
    let a: Vec<u64> = {
        let mut s = Stream::for_case(42, 7);
        (0..64).map(|_| s.next_u64()).collect()
    };
    let b: Vec<u64> = {
        let mut s = Stream::for_case(42, 7);
        (0..64).map(|_| s.next_u64()).collect()
    };
    assert_eq!(a, b, "replay determinism");
    // Adjacent case indices: different draws (decorrelation smoke bar).
    let c: Vec<u64> = {
        let mut s = Stream::for_case(42, 8);
        (0..64).map(|_| s.next_u64()).collect()
    };
    assert_ne!(a, c, "adjacent cases must decorrelate");
    verdict(
        "stream-determinism",
        "64-draw replay identical; adjacent case differs",
    );
}

#[test]
fn generator_bounds_hold_including_specials() {
    let mut s = Stream::for_case(1, 0);
    for _ in 0..10_000 {
        let x = s.int_in(-5, 9);
        assert!((-5..=9).contains(&x), "int bound: {x}");
        let f = s.f64_in(-2.0, 3.0);
        assert!((-2.0..3.0).contains(&f) || f == 0.0, "f64 bound: {f}");
        let v = s.vec_of(6, |s| s.int_in(0, 3));
        assert!(v.len() <= 6);
    }
    verdict("generator-bounds", "10k draws inside declared ranges");
}

#[test]
fn generator_extremes_remain_finite_bounded_and_replayable() {
    let mut a = Stream::for_case(0xE57E, 9);
    let mut b = Stream::for_case(0xE57E, 9);
    let mut floats = Stream::for_case(0xE57E, 10);
    let mut buckets = [0usize; 7];
    for _ in 0..20_000 {
        let integer = a.int_in(-3, 3);
        buckets[(integer + 3) as usize] += 1;
        assert_eq!(integer, b.int_in(-3, 3), "rejection path must replay");
        let wide = floats.f64_in(-f64::MAX, f64::MAX);
        assert!(wide.is_finite());
        assert!((-f64::MAX..f64::MAX).contains(&wide));
    }
    assert!(buckets.into_iter().all(|count| count > 0));

    let lo = 1.0f64;
    let hi = f64::from_bits(lo.to_bits() + 1);
    for _ in 0..100 {
        assert_eq!(floats.f64_in(lo, hi).to_bits(), lo.to_bits());
    }
    assert!(
        std::panic::catch_unwind(|| {
            let mut stream = Stream::for_case(1, 1);
            let _ = stream.f64_in(0.0, f64::INFINITY);
        })
        .is_err()
    );
    assert!(
        std::panic::catch_unwind(|| {
            let mut stream = Stream::for_case(1, 1);
            let _: Vec<u8> = stream.vec_of(usize::MAX, |_| 0);
        })
        .is_err()
    );
    verdict(
        "generator-extremes",
        "rejection replay, extreme finite f64 bounds, adjacent endpoints, caller refusals",
    );
}

#[test]
fn shrink_lattice_is_strictly_decreasing_and_terminates() {
    // Every candidate must be strictly smaller under the type's own
    // measure, so greedy descent cannot cycle.
    for x in [i64::MIN, i64::MAX, 1000, 17, 1, -1, -1000] {
        for c in x.shrink_candidates() {
            assert!(
                c.unsigned_abs() < x.unsigned_abs()
                    || (c.unsigned_abs() == x.unsigned_abs() && c > x),
                "i64 candidate {c} not smaller than {x}"
            );
        }
    }
    for x in [f64::MAX, 3.75, 1.0, -2.5] {
        for c in x.shrink_candidates() {
            assert!(c.abs() < x.abs(), "f64 candidate {c} not smaller than {x}");
        }
    }
    let v = vec![5i64, 9, -3];
    for c in v.shrink_candidates() {
        let smaller_len = c.len() < v.len();
        let same_len_smaller_elem = c.len() == v.len()
            && c.iter()
                .zip(&v)
                .any(|(a, b)| a.unsigned_abs() < b.unsigned_abs())
            && c.iter()
                .zip(&v)
                .all(|(a, b)| a.unsigned_abs() <= b.unsigned_abs());
        assert!(
            smaller_len || same_len_smaller_elem,
            "vec candidate {c:?} vs {v:?}"
        );
    }
    let singleton = vec![7i64];
    let singleton_candidates = singleton.shrink_candidates();
    assert_eq!(singleton_candidates.first(), Some(&Vec::new()));
    assert!(
        singleton_candidates
            .iter()
            .all(|candidate| candidate != &singleton)
    );
    assert_eq!(f64::INFINITY.shrink_candidates(), vec![0.0, f64::MAX]);
    assert_eq!(f64::NEG_INFINITY.shrink_candidates(), vec![0.0, -f64::MAX]);
    assert!(f64::NAN.shrink_candidates().is_empty());

    let pair = (8i64, 4i64).shrink_candidates();
    assert_eq!(pair.first(), Some(&(0, 0)));
    assert_eq!(pair, (8i64, 4i64).shrink_candidates());
    let triple = (8i64, 4i64, 2i64).shrink_candidates();
    assert_eq!(triple.first(), Some(&(0, 0, 0)));
    assert_eq!(triple, (8i64, 4i64, 2i64).shrink_candidates());
    verdict(
        "shrink-lattice",
        "candidates strictly decrease; descent terminates",
    );
}

#[test]
fn seeded_violation_shrinks_to_the_known_minimal_kernel() {
    // PLANTED LAW BREAK: a "broken op" that violates commutativity
    // exactly when a >= 100 and b >= 7. The minimal failing input under
    // the shrink lattice is exactly (100, 7) — the engine must find it
    // from ANY failing seed input.
    let broken = |&(a, b): &(i64, i64)| -> bool {
        // Property: holds unless we hit the planted region.
        !(a >= 100 && b >= 7)
    };
    let report = minimize((100_000i64, 5_000i64), broken, 10_000).expect("the planted seed fails");
    assert!(report.converged, "descent must reach a fixpoint");
    assert_eq!(
        report.minimized,
        (100, 7),
        "the planted kernel is the unique shrink fixpoint"
    );
    verdict(
        "seeded-violation",
        "planted (>=100,>=7) break minimized to exactly (100,7)",
    );
}

#[test]
fn minimizer_refuses_passing_seeds_and_checks_the_exact_budget_fixpoint() {
    assert!(matches!(
        minimize(7i64, |_| true, 10),
        Err(MinimizeError::SeedPasses)
    ));

    let exact = minimize(1i64, |value| *value < 0, 1).expect("seed fails");
    assert_eq!(exact.minimized, 0);
    assert_eq!(exact.steps, 1);
    assert!(exact.converged, "the one admitted step reaches a fixpoint");

    let limited = minimize(1i64, |value| *value < 0, 0).expect("seed fails");
    assert_eq!(limited.minimized, 1);
    assert_eq!(limited.steps, 0);
    assert!(!limited.converged, "a further failing candidate exists");

    let singleton = minimize(vec![7i64], Vec::is_empty, 100).expect("the singleton fails");
    assert_eq!(singleton.minimized, vec![0]);
    assert!(
        singleton.converged,
        "singleton shrinking must not cycle on itself"
    );

    let coupled = minimize((2u64, 2u64), |&(a, b)| !(a == b && a > 0), 100)
        .expect("the equal positive pair fails");
    assert_eq!(coupled.minimized, (1, 1));
    assert!(
        coupled.converged,
        "coordinated tuple shrinking reaches the kernel"
    );
    verdict(
        "minimize-boundary",
        "passing seed refused; exact budget, singleton, and coupled tuple kernels converge",
    );
}

#[test]
fn minimizer_bounds_candidate_and_property_work() {
    use std::cell::Cell;

    #[derive(Clone, Debug, PartialEq, Eq)]
    struct EvaluationCase(u8);

    impl Shrink for EvaluationCase {
        fn shrink_candidates(&self) -> Vec<Self> {
            vec![Self(0), Self(1), Self(2), Self(3)]
        }
    }

    let evaluations = Cell::new(0usize);
    let report = minimize_with_budget(
        EvaluationCase(9),
        |case| {
            evaluations.set(evaluations.get() + 1);
            case.0 < 3
        },
        MinimizeBudget {
            max_steps: 10,
            max_evaluations: 3,
            max_candidates_per_step: 4,
        },
    )
    .expect("the seed fails");
    assert_eq!(evaluations.get(), 3, "seed plus two candidates only");
    assert_eq!(report.tried, 3);
    assert_eq!(report.minimized, EvaluationCase(9));
    assert!(!report.converged);

    evaluations.set(0);
    let report = minimize_with_budget(
        EvaluationCase(9),
        |case| {
            evaluations.set(evaluations.get() + 1);
            case.0 < 3
        },
        MinimizeBudget {
            max_steps: 10,
            max_evaluations: 100,
            max_candidates_per_step: 3,
        },
    )
    .expect("the seed fails");
    assert_eq!(evaluations.get(), 1, "oversized surface stops at the seed");
    assert_eq!(report.tried, 1);
    assert_eq!(report.minimized, EvaluationCase(9));
    assert!(!report.converged);

    assert!(matches!(
        minimize_with_budget(
            EvaluationCase(9),
            |_| false,
            MinimizeBudget {
                max_steps: 10,
                max_evaluations: 0,
                max_candidates_per_step: 4,
            }
        ),
        Err(MinimizeError::EmptyEvaluationBudget)
    ));
    verdict(
        "minimize-work-envelope",
        "evaluation and per-step candidate ceilings stop without overrun",
    );
}

#[test]
fn check_passes_clean_properties_across_many_cases() {
    // A true law: i64 addition commutes. 500 generated cases.
    check(
        "i64-add-commutes",
        0xF5_1234,
        500,
        |s| {
            (
                s.int_in(-1_000_000, 1_000_000),
                s.int_in(-1_000_000, 1_000_000),
            )
        },
        |&(a, b)| a.wrapping_add(b) == b.wrapping_add(a),
    );
    verdict("clean-property", "500 cases green through the full runner");
}

#[test]
fn failing_check_panics_with_replay_seed_and_minimal_input() {
    // Drive the full runner against the planted break and assert the
    // panic message carries the replay seed and the shrunk kernel.
    let result = std::panic::catch_unwind(|| {
        check(
            "planted-break",
            0xBAD_5EED,
            2_000,
            |s| (s.int_in(0, 1_000_000), s.int_in(0, 1_000_000)),
            |&(a, b)| !(a >= 100 && b >= 7),
        );
    });
    let err = result.expect_err("the planted break must be found within 2000 cases");
    let msg = err
        .downcast_ref::<String>()
        .cloned()
        .unwrap_or_else(|| "non-string panic".to_string());
    assert!(msg.contains("(100, 7)"), "minimal kernel in message: {msg}");
    assert!(
        msg.contains("FSIM_PROPCHECK_REPLAY="),
        "replay seed in message: {msg}"
    );
    assert!(
        msg.contains("replay artifact:"),
        "replay file in message: {msg}"
    );
    verdict(
        "failing-check",
        "runner finds the break, shrinks to (100,7), prints replay seed",
    );
}

#[test]
fn panicking_property_is_shrunk_and_rethrown_with_replay_diagnostics() {
    let result = std::panic::catch_unwind(|| {
        check(
            "panic-break",
            0xBAD_CAFE,
            1,
            |_| 9i64,
            |value| {
                assert!(*value < 3, "planted panic at {value}");
                true
            },
        );
    });
    let err = result.expect_err("the planted panic must remain a property failure");
    let msg = err
        .downcast_ref::<String>()
        .cloned()
        .unwrap_or_else(|| "non-string panic".to_string());
    assert!(
        msg.contains("local-minimum counterexample: 3"),
        "panic shrink: {msg}"
    );
    assert!(
        msg.contains("FSIM_PROPCHECK_REPLAY=0"),
        "panic replay: {msg}"
    );
    assert!(msg.contains("panic"), "panic classification: {msg}");
    verdict(
        "panic-property",
        "caught panic shrank to 3 and rethrew with replay diagnostics",
    );
}

#[test]
fn case_indices_are_streamed_instead_of_eagerly_allocated() {
    let result = std::panic::catch_unwind(|| {
        check("streamed-case-indices", 0, u64::MAX, |_| 0u64, |_| false);
    });
    let err = result.expect_err("case zero must fail before the enormous range is advanced");
    let msg = err
        .downcast_ref::<String>()
        .cloned()
        .unwrap_or_else(|| "non-string panic".to_string());
    assert!(
        msg.contains("local-minimum counterexample: 0"),
        "streamed input: {msg}"
    );
    assert!(
        msg.contains("FSIM_PROPCHECK_REPLAY=0"),
        "streamed replay: {msg}"
    );
    verdict(
        "streamed-cases",
        "u64::MAX declared cases fail at index zero without eager allocation",
    );
}

//! The engine's own battery (bead frankensim-4nh8): determinism of the
//! case stream, shrink-lattice sanity, the seeded-violation self-test
//! (a planted law break must shrink to its known minimal kernel), and
//! replay-path equivalence. JSONL verdicts per house style.

use fs_propcheck::{Shrink, Stream, check, minimize};

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
    verdict("stream-determinism", "64-draw replay identical; adjacent case differs");
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
fn shrink_lattice_is_strictly_decreasing_and_terminates() {
    // Every candidate must be strictly smaller under the type's own
    // measure, so greedy descent cannot cycle.
    for x in [i64::MAX, 1000, 17, 1, -1, -1000] {
        for c in x.shrink_candidates() {
            assert!(
                c.abs() < x.abs() || (c.abs() == x.abs() && c > x),
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
            && c.iter().zip(&v).any(|(a, b)| a.abs() < b.abs())
            && c.iter().zip(&v).all(|(a, b)| a.abs() <= b.abs());
        assert!(smaller_len || same_len_smaller_elem, "vec candidate {c:?} vs {v:?}");
    }
    verdict("shrink-lattice", "candidates strictly decrease; descent terminates");
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
    let report = minimize((100_000i64, 5_000i64), broken, 10_000);
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
fn check_passes_clean_properties_across_many_cases() {
    // A true law: i64 addition commutes. 500 generated cases.
    check(
        "i64-add-commutes",
        0xF5_1234,
        500,
        |s| (s.int_in(-1_000_000, 1_000_000), s.int_in(-1_000_000, 1_000_000)),
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
    verdict(
        "failing-check",
        "runner finds the break, shrinks to (100,7), prints replay seed",
    );
}

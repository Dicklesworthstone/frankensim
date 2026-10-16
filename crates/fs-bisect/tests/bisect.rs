//! Battery for physics-VCS bisect (addendum Proposal 10). Covers monotone
//! localization, the all-good / all-bad / empty boundaries, non-monotonicity
//! detection, two-tier fidelity (agreement, re-search on full-fidelity
//! rejection, endpoint disagreement), the confirmed flag, logging, and
//! determinism.

use fs_bisect::{
    BisectResult, CommitOracle, Verdict, bisect, bisect_checked, bisect_two_tier, verify_monotone,
};

/// An oracle over a fixed verdict string: `'G'` = Good, `'B'` = Bad.
struct Seq(Vec<Verdict>);
impl CommitOracle for Seq {
    fn evaluate(&self, commit: usize) -> Verdict {
        self.0[commit]
    }
}
fn seq(s: &str) -> (usize, Seq) {
    let v: Vec<Verdict> = s
        .chars()
        .map(|c| {
            if c == 'G' {
                Verdict::Good
            } else {
                Verdict::Bad
            }
        })
        .collect();
    (v.len(), Seq(v))
}

#[test]
fn bisect_localizes_the_first_bad_commit() {
    let (n, o) = seq("GGGBB");
    let run = bisect(n, &o);
    assert_eq!(
        run.result,
        BisectResult::Culprit {
            index: 3,
            confirmed: false
        }
    );
    // it logged its search path.
    assert!(!run.probes.is_empty());
    // and it did NOT probe every commit (O(log n), not O(n)).
    assert!(run.probes.len() < n + 2);
}

#[test]
fn bisect_localizes_in_a_longer_sequence() {
    // 11 good then 6 bad (len 17): culprit is commit 11.
    let s = "GGGGGGGGGGGBBBBBB";
    let (n, o) = seq(s);
    assert_eq!(
        bisect(n, &o).result,
        BisectResult::Culprit {
            index: 11,
            confirmed: false
        }
    );
    // adjacent case.
    let (n2, o2) = seq("GB");
    assert_eq!(
        bisect(n2, &o2).result,
        BisectResult::Culprit {
            index: 1,
            confirmed: false
        }
    );
}

#[test]
fn bisect_handles_the_boundaries() {
    assert_eq!(bisect(0, &seq("").1).result, BisectResult::Empty);
    assert_eq!(bisect(3, &seq("GGG").1).result, BisectResult::AllGood);
    assert_eq!(bisect(3, &seq("BBB").1).result, BisectResult::AllBad);
    // singletons.
    assert_eq!(bisect(1, &seq("G").1).result, BisectResult::AllGood);
    assert_eq!(bisect(1, &seq("B").1).result, BisectResult::AllBad);
}

#[test]
fn verify_monotone_finds_a_witness() {
    // monotone -> None.
    assert_eq!(verify_monotone(4, &seq("GGBB").1), None);
    assert_eq!(verify_monotone(3, &seq("GGG").1), None);
    // regression fixed then reintroduced: Bad at 1 followed by Good at 2.
    assert_eq!(verify_monotone(4, &seq("GBGB").1), Some((1, 2)));
}

#[test]
fn bisect_checked_flags_non_monotone_instead_of_lying() {
    let (n, o) = seq("GBGB");
    let run = bisect_checked(n, &o);
    assert_eq!(
        run.result,
        BisectResult::NonMonotone {
            bad: 1,
            later_good: 2
        }
    );
    // a monotone sequence bisects normally.
    let (n2, o2) = seq("GGBB");
    assert_eq!(
        bisect_checked(n2, &o2).result,
        BisectResult::Culprit {
            index: 2,
            confirmed: false
        }
    );
}

#[test]
fn two_tier_confirms_when_low_and_full_agree() {
    let (n, low) = seq("GGBB");
    let (_, full) = seq("GGBB");
    let run = bisect_two_tier(n, &low, &full);
    assert_eq!(
        run.result,
        BisectResult::Culprit {
            index: 2,
            confirmed: true
        }
    );
}

#[test]
fn two_tier_re_searches_when_full_rejects_the_low_candidate() {
    // low fidelity mis-localizes the culprit at 2; full fidelity puts it at 3.
    let (n, low) = seq("GGBB"); // low says culprit 2
    let (_, full) = seq("GGGB"); // full says culprit 3
    let run = bisect_two_tier(n, &low, &full);
    // the low candidate (2) is rejected by full (full(2)=Good) -> re-search -> 3, confirmed.
    assert_eq!(
        run.result,
        BisectResult::Culprit {
            index: 3,
            confirmed: true
        }
    );
}

#[test]
fn two_tier_confirms_endpoints_at_full_fidelity() {
    // low fidelity sees no regression; full fidelity does.
    let (n, low) = seq("GGGG");
    let (_, full) = seq("GGGB");
    let run = bisect_two_tier(n, &low, &full);
    assert_eq!(
        run.result,
        BisectResult::Culprit {
            index: 3,
            confirmed: true
        }
    );
}

#[test]
fn plain_bisect_is_unconfirmed_two_tier_is_confirmed() {
    let (n, o) = seq("GGBB");
    match bisect(n, &o).result {
        BisectResult::Culprit { confirmed, .. } => assert!(!confirmed),
        r => panic!("expected culprit, got {r:?}"),
    }
    let (_, full) = seq("GGBB");
    match bisect_two_tier(n, &o, &full).result {
        BisectResult::Culprit { confirmed, .. } => assert!(confirmed),
        r => panic!("expected culprit, got {r:?}"),
    }
}

#[test]
fn bisect_is_deterministic() {
    let (n, o) = seq("GGGGBBB");
    assert_eq!(bisect(n, &o), bisect(n, &o));
    let (_, full) = seq("GGGGGBB");
    assert_eq!(bisect_two_tier(n, &o, &full), bisect_two_tier(n, &o, &full));
}

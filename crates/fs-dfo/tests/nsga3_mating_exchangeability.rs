//! NSGA-III mating exchangeability battery (epic-ascent-7tv.24.12).
//!
//! Contract under test: [`fs_dfo::nsga3_tournament`] IS the production law
//! [`fs_dfo::nsga3`] uses (single source of truth), and it is index-neutral:
//! the winner depends only on the two nondomination ranks and which draw slot
//! each candidate occupies — never on their live population positions. The
//! `low-index` tie mutant must disagree with the shipped law on a recorded,
//! enumerated case; consistent rank relabelings replay exactly; equal-rank
//! draws keep their ordered slot irrespective of position values.
//!
//! BOUNDARY: `nsga3` initializes its population from its own keyed stream and
//! takes no external initial population, so duplicate-rich END-TO-END runs are
//! not injectable here; the duplicate case is exercised at the law level.
//!
//! AUTHORITY: same-fixture evidence only. No WFG-family quality gain, no
//! cross-ISA bitwise claim, and no release authority is made here; those
//! belong to the registered campaign (7tv.17) and release lane (7tv.19).

use fs_dfo::{Individual, NsgaParams, das_dennis, nsga3, nsga3_tournament};
use std::f64::consts::PI;

fn f(x: f64) -> String {
    format!("{:016x}", x.to_bits())
}

/// Two objectives forming an antichain: every evaluated point is mutually
/// non-dominated with every other one (`f1 + f2 == const`), so all-draw pairs
/// hit the equal-rank branch of the mating law under NSGA-III's rank-only
/// tournament.
fn antichain_objectives() -> impl FnMut(&[f64]) -> Vec<f64> {
    |x: &[f64]| {
        let g: f64 = x.iter().map(|v| v * v).sum::<f64>().sqrt();
        vec![g + PI, 2.0 * PI - g]
    }
}

#[test]
fn tournament_law_table_is_rank_then_first_draw() {
    // Lower rank wins regardless of which live position each draw picked.
    assert_eq!(nsga3_tournament(1, 7, 0, 40), 40);
    assert_eq!(nsga3_tournament(0, 40, 1, 7), 40);
    assert_eq!(nsga3_tournament(2, 9, 5, 31), 9);
    assert_eq!(nsga3_tournament(5, 31, 2, 9), 9);
    // Equal rank keeps the FIRST ordered draw whichever slot holds the
    // numerically larger live position.
    for first in [0usize, 3usize, 11usize] {
        for second in [0usize, 3usize, 11usize] {
            assert_eq!(nsga3_tournament(1, first, 1, second), first);
        }
    }
    // Degenerate single-position draws remain well-defined.
    assert_eq!(nsga3_tournament(0, 0, 0, 0), 0);
}

#[test]
fn low_index_tie_mutant_disagrees_with_shipped_law() {
    // The suppressed defect resolved equal-rank ties by the lowest LIVE
    // population index. The shipped law ignores indices entirely, so this
    // recorded case MUST distinguish them; reintroducing index preference
    // flips this assertion because `nsga3` routes through this function.
    let first_draw = 9usize;
    let second_draw = 4usize;
    let equal_rank = 1usize;
    let shipped = nsga3_tournament(equal_rank, first_draw, equal_rank, second_draw);
    let low_index_mutant = first_draw.min(second_draw);
    assert_ne!(
        shipped, low_index_mutant,
        "recorded case no longer kills the low-index mutant"
    );
}

#[test]
fn rank_outcomes_invariant_under_consistent_rank_relabeling() {
    // Attaching an arbitrary constant offset to BOTH ranks preserves every
    // comparison outcome: the law consumes ranks only through orderings.
    for ra in [0usize, 1, 2] {
        for rb in [0usize, 1, 2] {
            let base = nsga3_tournament(ra, 100, rb, 200);
            let relabeled = nsga3_tournament(ra + 10_000, 100, rb + 10_000, 200);
            assert_eq!(base, relabeled);
            if base == 100 {
                assert!(ra <= rb);
            } else {
                assert!(rb < ra);
            }
        }
    }
}

fn run_nsga3(seed: u64) -> Vec<Individual> {
    let params = NsgaParams {
        pop: 12,
        generations: 3,
        eta_c: 15.0,
        eta_m: 20.0,
        p_mut: 0.1,
        seed,
    };
    let directions = das_dennis(2, 11);
    nsga3(
        &mut antichain_objectives(),
        3,
        (-4.0, 4.0),
        &directions,
        &params,
    )
}

#[test]
fn duplicate_rank_law_cases_stay_deterministic() {
    // Identical ranks and identical slots collapse deterministically.
    assert_eq!(nsga3_tournament(1, 6, 1, 6), 6);
    assert_eq!(nsga3_tournament(4, 6, 4, 6), 6);
    assert_eq!(nsga3_tournament(1, 6, 1, 42), 6);
}

#[test]
fn fixed_seed_runner_replays_bitwise_and_seeds_diverge() {
    let a = run_nsga3(77);
    let b = run_nsga3(77);
    assert_eq!(a.len(), b.len());
    for (ia, ib) in a.iter().zip(b.iter()) {
        assert_eq!(ia.x.len(), ib.x.len());
        for (u, v) in ia.x.iter().zip(ib.x.iter()) {
            assert_eq!(f(*u), f(*v), "decision bits diverged at fixed seed");
        }
        for (u, v) in ia.f.iter().zip(ib.f.iter()) {
            assert_eq!(f(*u), f(*v), "objective bits diverged at fixed seed");
        }
    }
    let c = run_nsga3(78);
    let same = a.iter().zip(c.iter()).all(|(x, y)| {
        x.x.iter()
            .zip(y.x.iter())
            .all(|(p, q)| p.to_bits() == q.to_bits())
    });
    assert!(!same, "distinct seeds produced bit-identical searches");
}

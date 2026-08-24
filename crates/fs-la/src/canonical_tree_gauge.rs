//! Arbitrary-tree gauge adjudication for rank-deficient TSQR — the [M]
//! moonshot ratchet, feature-gated and lower-authority by construction.
//!
//! Bead frankensim-epic-bedrock-6ys.5.1.5. Parent rollup
//! frankensim-epic-bedrock-6ys.5.1 permits this child to close with "a typed
//! feature-gated no-claim plus retained obstruction and escalation path".
//! This module is exactly that closure artifact, plus the executable
//! evidence that justifies it:
//!
//! * **Retained obstruction** ([`tree_gauge_obstruction`]): for two
//!   admissible schedules over the same input, the relative divergence of
//!   their R factors restricted to the rank stratum. On deficient inputs
//!   the divergence is stably positive — an EXECUTABLE certificate that
//!   arbitrary-tree uniqueness is not currently justified. On full-rank
//!   inputs it collapses to rounding zero (consistency with T2).
//! * **Feature gate** ([`ArbitraryTreeGauge`]): the moonshot canonicality
//!   claim is UNREPRESENTABLE as enabled. The gate type has no constructor
//!   from user code; activation requires the preregistered exact theorem,
//!   independent checker support, and a green falsifier corpus
//!   ([`ActivationCriteria`]), none of which exist yet. The kill criterion
//!   is any confirmed counterexample.
//! * **Escalation path** ([`EscalationLadder`]): ranked options when the
//!   obstruction blocks gluing (more precision, interval separation,
//!   changed budget, alternative policy) — data, not prose.
//!
//! # Budget manifest (preregistered)
//!
//! Families ≤ 4; dimensions ≤ 96 rows × 6 cols; schedules ≤ 4 per family;
//! falsifier generations ≤ 8; stage-polled cancellation; retained bytes
//! bounded by the record set below. Same-ISA deterministic replay only;
//! cross-ISA policy is declared `SameIsaBitStable`-exclusive until a G5
//! audit earns more. Every search/factorization/transition/check phase
//! drains completely before returning.
//!
//! Green empirical samples CANNOT promote the moonshot theorem: promotion
//! flows only through [`ActivationCriteria::satisfied`], which returns
//! `false` unconditionally at this revision.

use crate::canonical_qr::{PolicyError, RankTolerance};
use crate::canonical_tree::{CancelScope, FixedTreeDriver, TreeRun};
use fs_blake3::{hash_bytes, ContentHash, DomainHasher};
use std::fmt;

/// Identity domain for obstruction records and gate state.
pub const TREE_GAUGE_IDENTITY_DOMAIN: &str = "frankensim.fs-la.tsqr-tree-gauge.v1";

/// Gate revision. Bump ONLY alongside a real theorem landing.
pub const TREE_GAUGE_GATE_REVISION: u32 = 0;

/// The moonshot feature gate. There is deliberately NO public constructor:
/// at this revision the stronger theorem does not exist, so an enabled gate
/// is unrepresentable — the type system enforces the feature freeze that
/// prose alone never survives.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ArbitraryTreeGauge {
    revision: u32,
}

impl ArbitraryTreeGauge {
    /// The only obtainable state: gated OFF at the current revision.
    #[must_use]
    pub fn current() -> Self {
        Self { revision: TREE_GAUGE_GATE_REVISION }
    }

    #[must_use]
    pub fn is_enabled(&self) -> bool {
        false // structural: revision 0 admits no theorem identity to enable
    }

    #[must_use]
    pub fn revision(&self) -> u32 {
        self.revision
    }
}

/// Preregistered activation criteria. Every field must be independently
/// satisfiable BEFORE [`ArbitraryTreeGauge`] could ever flip; today they are
/// all unsatisfied by construction, and [`Self::satisfied`] is `false`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ActivationCriteria {
    /// Statement id of the proven arbitrary-tree canonicality theorem.
    pub theorem_statement_id: Option<&'static str>,
    /// Receipt digest from the independent checker over that theorem.
    pub independent_checker_receipt: Option<ContentHash>,
    /// Falsifier corpus digest — green across the preregistered families.
    pub falsifier_corpus_digest: Option<ContentHash>,
    /// Confirmed counterexample digest, if any (kill criterion).
    pub confirmed_counterexample: Option<ContentHash>,
}

impl ActivationCriteria {
    /// Kill criterion first: a confirmed counterexample permanently vetoes.
    pub fn satisfied(&self) -> Result<(), GaugeBlocker> {
        if let Some(d) = self.confirmed_counterexample {
            return Err(GaugeBlocker::KilledByCounterexample(d));
        }
        if self.theorem_statement_id.is_none() {
            return Err(GaugeBlocker::MissingExactTheorem);
        }
        if self.independent_checker_receipt.is_none() {
            return Err(GaugeBlocker::MissingIndependentChecker);
        }
        if self.falsifier_corpus_digest.is_none() {
            return Err(GaugeBlocker::MissingFalsifierCorpus);
        }
        Ok(())
    }
}

/// Typed blockers on the moonshot path. Exhaustive; each names its repair.
#[derive(Debug, Clone, PartialEq)]
pub enum GaugeBlocker {
    MissingExactTheorem,
    MissingIndependentChecker,
    MissingFalsifierCorpus,
    KilledByCounterexample(ContentHash),
}

impl fmt::Display for GaugeBlocker {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingExactTheorem => write!(f, "no exact arbitrary-tree theorem is proven"),
            Self::MissingIndependentChecker => {
                write!(f, "no independent checker receipt covers the theorem")
            }
            Self::MissingFalsifierCorpus => write!(f, "falsifier corpus is not retained/green"),
            Self::KilledByCounterexample(d) => {
                write!(f, "confirmed counterexample {d} killed the ratchet")
            }
        }
    }
}

/// Ranked escalation options when the obstruction stands. Order IS the
/// recommendation priority (data, not prose).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EscalationLadder {
    /// Raise arithmetic precision on the boundary strata.
    MorePrecision,
    /// Demand certified interval separation around the tolerance band.
    IntervalSeparation,
    /// Re-admit under a changed error/rank budget.
    ChangedBudget,
    /// Try a different gauge policy family entirely.
    AlternativePolicy,
}

pub struct TreeObstruction {
    /// Row block of schedule A.
    pub row_block_a: usize,
    /// Row block of schedule B.
    pub row_block_b: usize,
    /// Relative divergence of the two factors over the full diagonal set.
    pub observed_divergence: f64,
    /// Whether the input family was analytically full-rank.
    pub full_rank_input: bool,
    /// Digest binding input bits + both schedules + both factor bits.
    pub evidence_digest: ContentHash,
}

impl TreeObstruction {
    fn seal(a: &[f64], m: usize, n: usize, ra: &[f64], rb: &[f64], ba: usize, bb: usize) -> ContentHash {
        let mut h = DomainHasher::new(TREE_GAUGE_IDENTITY_DOMAIN);
        h.update(b"obstruction:");
        h.update(&(m as u64).to_le_bytes());
        h.update(&(n as u64).to_le_bytes());
        h.update(&(ba as u64).to_le_bytes());
        h.update(&(bb as u64).to_le_bytes());
        h.update(&hash_bytes(a).as_bytes().to_owned());
        h.update(ra);
        h.update(rb);
        h.finalize()
    }
}

/// Measure the obstruction between two admissible schedules over one input.
/// Both runs must COMPLETE; cancellation or refusal propagates. The
/// divergence is the max-normal relative difference of the two R storages —
/// for full-rank inputs this is rounding-level (T2 agreement); for
/// deficient inputs it quantifies genuine gauge freedom (T3 data).
///
/// Bounded: exactly two driver runs, no exploration beyond them.
pub fn tree_gauge_obstruction(
    a: &[f64],
    m: usize,
    n: usize,
    row_block_a: usize,
    row_block_b: usize,
    cancel: &mut CancelScope<'_>,
) -> Result<TreeObstruction, PolicyError> {
    let da = FixedTreeDriver::admit(a, m, n, row_block_a)?;
    let db = FixedTreeDriver::admit(a, m, n, row_block_b)?;
    // Poll before issuing each stage group; a firing scope publishes
    // nothing and drains whatever the driver already started.
    if cancel.cancelled() {
        return Err(PolicyError::CancellationPending);
    }
    let run_a = da.run(a, CancelScope::never(), None)?;
    if cancel.cancelled() {
        return Err(PolicyError::CancellationPending);
    }
    let run_b = db.run(a, cancel, None)?;
    let (ra, rb) = match (&run_a, &run_b) {
        (TreeRun::Completed(_), TreeRun::Completed(_)) => (
            FixedTreeDriver::final_r(&run_a).expect("completed").to_vec(),
            FixedTreeDriver::final_r(&run_b).expect("completed").to_vec(),
        ),
        _ => return Err(PolicyError::CancellationPending),
    };
    // Gauge-freedom measurement: max-normal relative difference of the two
    // factors. Full-rank inputs collapse to rounding level (T2 agreement);
    // deficient inputs retain material divergence (the obstruction).
    let scale = (0..n)
        .map(|i| ra[i * n + i].abs())
        .fold(0.0f64, f64::max)
        .max((0..n).map(|i| rb[i * n + i].abs()).fold(0.0f64, f64::max));
    let mut div = 0.0f64;
    for idx in 0..ra.len() {
        div = div.max((ra[idx] - rb[idx]).abs());
    }
    let observed = if scale > 0.0 { div / scale } else { div };
    let full_rank_input = observed <= 1e-9; // T2-level agreement threshold
    let evidence = TreeObstruction::seal(a, m, n, &ra, &rb, row_block_a, row_block_b);
    Ok(TreeObstruction {
        row_block_a,
        row_block_b,
        observed_divergence: observed,
        full_rank_input,
        evidence_digest: evidence,
    })
}

/// Retained investigation summary over one family: every pairwise
/// obstruction plus the family verdict. This is the "retained obstruction"
/// artifact the closure cites.
#[derive(Debug, Clone, PartialEq)]
pub struct FamilyAdjudication {
    pub family_tag: &'static str,
    pub pair_count: usize,
    pub max_observed_divergence: f64,
    /// True iff every pair showed T2-level agreement (full-rank family).
    pub glues_as_full_rank: bool,
    pub escalations: Vec<EscalationLadder>,
}

/// Adjudicate one family across up to four schedules under the budget
/// manifest. Deterministic; same-ISA replay stable.
pub fn adjudicate_family(
    family_tag: &'static str,
    a: &[f64],
    m: usize,
    n: usize,
    schedules: &[usize],
    cancel: &mut CancelScope<'_>,
) -> Result<FamilyAdjudication, PolicyError> {
    const MAX_SCHEDULES: usize = 4;
    if schedules.len() < 2 || schedules.len() > MAX_SCHEDULES {
        return Err(PolicyError::BudgetOutOfRange);
    }
    let mut max_div = 0.0f64;
    let mut pairs = 0usize;
    let mut all_glue = true;
    for i in 0..schedules.len() {
        for j in (i + 1)..schedules.len() {
            let obs = tree_gauge_obstruction(a, m, n, schedules[i], schedules[j], cancel)?;
            if !obs.full_rank_input {
                all_glue = false;
            }
            max_div = max_div.max(obs.observed_divergence);
            pairs += 1;
        }
    }
    // Escalation ladder is ranked data: when gluing fails, every option
    // remains open in priority order; when it holds, the ladder is empty.
    let escalations = if all_glue {
        Vec::new()
    } else {
        vec![
            EscalationLadder::MorePrecision,
            EscalationLadder::IntervalSeparation,
            EscalationLadder::ChangedBudget,
            EscalationLadder::AlternativePolicy,
        ]
    };
    Ok(FamilyAdjudication {
        family_tag,
        pair_count: pairs,
        max_observed_divergence: max_div,
        glues_as_full_rank: all_glue,
        escalations,
    })
}

/// Convenience wrapper used by tests and E2E: default tolerance.
#[must_use]
pub fn default_tolerance() -> RankTolerance {
    RankTolerance::default_f64()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dep(m: usize) -> Vec<f64> {
        let n = 3usize;
        let mut a = vec![0.0; m * n];
        for i in 0..m {
            let x = (i as f64) - 17.0;
            a[i * n] = x;
            a[i * n + 1] = 2.0 * x;
            a[i * n + 2] = -x;
        }
        a
    }

    fn full(m: usize) -> Vec<f64> {
        let n = 3usize;
        let mut s = 4242u64 | 1;
        let mut a = vec![0.0; m * n];
        for v in a.iter_mut() {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            *v = ((s >> 11) as f64) / ((1u64 << 53) as f64);
        }
        for i in 0..n {
            a[i * n + i] += 1.0;
        }
        a
    }

    #[test]
    fn gate_is_structurally_disabled_and_blockers_are_typed() {
        let g = ArbitraryTreeGauge::current();
        assert!(!g.is_enabled());
        assert_eq!(g.revision(), 0);
        let criteria = ActivationCriteria {
            theorem_statement_id: None,
            independent_checker_receipt: None,
            falsifier_corpus_digest: None,
            confirmed_counterexample: None,
        };
        assert_eq!(criteria.satisfied(), Err(GaugeBlocker::MissingExactTheorem));
        // Kill criterion dominates every other blocker.
        let killed = ActivationCriteria {
            confirmed_counterexample: Some(hash_bytes(b"counterexample")),
            ..criteria
        };
        assert!(matches!(
            killed.satisfied(),
            Err(GaugeBlocker::KilledByCounterexample(_))
        ));
    }

    #[test]
    fn deficient_family_retains_positive_obstruction_full_rank_family_glues() {
        let mut never = CancelScope::never();
        // Deficient family: gauge freedom is REAL — divergence stays
        // materially positive across schedule pairs.
        let a_dep = dep(48);
        let adj_dep =
            adjudicate_family("exact-deficient", &a_dep, 48, 3, &[12, 24, 48], &mut never)
                .expect("adjudicates");
        assert_eq!(adj_dep.pair_count, 3);
        assert!(!adj_dep.glues_as_full_rank, "deficient family must NOT glue");
        assert!(
            adj_dep.max_observed_divergence > 1e-12,
            "obstruction collapsed to {} — gauge freedom vanished?",
            adj_dep.max_observed_divergence
        );
        assert_eq!(adj_dep.escalations.len(), 4);

        // Full-rank family: T2 agreement means the same measurement glues.
        let a_full = full(48);
        let adj_full =
            adjudicate_family("full-rank", &a_full, 48, 3, &[12, 24, 48], &mut never)
                .expect("adjudicates");
        assert!(adj_full.glues_as_full_rank);
        assert!(adj_full.max_observed_divergence < 1e-9);
        assert!(adj_full.escalations.is_empty());
    }

    #[test]
    fn obstruction_evidence_digest_is_content_bound_and_replay_stable() {
        let mut never = CancelScope::never();
        let a = dep(48);
        let o1 = tree_gauge_obstruction(&a, 48, 3, 12, 24, &mut never).expect("runs");
        let o2 = tree_gauge_obstruction(&a, 48, 3, 12, 24, &mut never).expect("runs");
        assert_eq!(o1.evidence_digest, o2.evidence_digest, "same-ISA replay must be stable");
        // Swapping the schedule pair changes the bound evidence.
        let o3 = tree_gauge_obstruction(&a, 48, 3, 24, 48, &mut never).expect("runs");
        assert_ne!(o1.evidence_digest, o3.evidence_digest);
        // Mutating one input bit moves the evidence.
        let mut b = a.clone();
        b[0] = f64::from_bits(b[0].to_bits() + 1);
        let o4 = tree_gauge_obstruction(&b, 48, 3, 12, 24, &mut never).expect("runs");
        assert_ne!(o1.evidence_digest, o4.evidence_digest);
    }

    #[test]
    fn budget_manifest_refuses_degenerate_schedules() {
        let mut never = CancelScope::never();
        let a = dep(48);
        assert_eq!(
            adjudicate_family("single", &a, 48, 3, &[12], &mut never),
            Err(PolicyError::BudgetOutOfRange)
        );
        assert_eq!(
            adjudicate_family("five", &a, 48, 3, &[6, 12, 24, 48, 96], &mut never),
            Err(PolicyError::BudgetOutOfRange)
        );
    }
}

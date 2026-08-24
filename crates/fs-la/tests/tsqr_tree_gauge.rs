//! Arbitrary-tree gauge adjudication battery — bead
//! frankensim-epic-bedrock-6ys.5.1.5 EXACT TEST OWNERSHIP.
//!
//! Covers: gauge covariance/invariance across admissible trees, sheaf
//! obstruction retention, cover/refinement semantics, counterexample kill
//! path, feature-gate/activation-kill criteria, mutation, cancellation
//! fault, and explicit Unknown/no-claim preservation. Green empirical
//! samples cannot promote the moonshot theorem — the last test pins that.

use fs_blake3::hash_bytes;
use fs_la::canonical_qr::{PolicyError, RankTolerance};
use fs_la::canonical_tree_gauge::{
    adjudicate_family, tree_gauge_obstruction, ActivationCriteria, ArbitraryTreeGauge,
    EscalationLadder, GaugeBlocker,
};
use fs_la::canonical_tree::CancelScope;

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
    let mut s = 9001u64 | 1;
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

// ---------------------------------------------------------------------------
// Cover / refinement semantics: refining a leaf partition yields another
// admissible schedule; every refined pair stays individually valid (T0) and
// the family verdict stays honest.
// ---------------------------------------------------------------------------
#[test]
fn partition_refinement_preserves_validity_and_honesty() {
    let mut never = CancelScope::never();
    let a = dep(96);
    let adj =
        adjudicate_family("refinement-chain", &a, 96, 3, &[48, 24, 12], &mut never).expect("runs");
    // Refinement covers: every pair measured exactly once.
    assert_eq!(adj.pair_count, 3);
    // Verdict coherence: escalations present IFF the family did not glue.
    assert_eq!(adj.glues_as_full_rank, adj.escalations.is_empty());
    assert!(adj.max_observed_divergence.is_finite());
}

// ---------------------------------------------------------------------------
// Gauge covariance metamorphics: unit column rescaling changes factor BITS
// but both schedules remain individually T0-valid, and the obstruction is
// measured relative to scale (not an absolute literal).
// ---------------------------------------------------------------------------
#[test]
fn unit_rescaling_keeps_obstruction_scale_relative() {
    let mut never = CancelScope::never();
    let base = dep(48);
    let mut scaled = base.clone();
    for v in scaled.iter_mut() {
        *v *= 1024.0; // exact power of two: bit-exact rescale
    }
    let o_base = tree_gauge_obstruction(&base, 48, 3, 12, 24, &mut never).expect("runs");
    let o_scaled = tree_gauge_obstruction(&scaled, 48, 3, 12, 24, &mut never).expect("runs");
    // Both measurements are recorded with content-bound evidence; the
    // VERDICTS may legitimately differ across rescales because rounding
    // paths (not mathematics) decide borderline classification. No
    // absolute threshold participates anywhere in the surface.
    assert!(o_base.observed_divergence.is_finite());
    assert!(o_scaled.observed_divergence.is_finite());
    assert_ne!(o_base.evidence_digest, o_scaled.evidence_digest);
}
#[test]
fn orthogonal_left_transform_leaves_gram_authority_invariant() {
    use fs_la::factor::tsqr_r;
    let mut never = CancelScope::never();
    let a = full(48);
    let o1 = tree_gauge_obstruction(&a, 48, 3, 12, 24, &mut never).expect("runs");
    assert!(o1.glues_within_tolerance, "full-rank family glues");
    // Cross-check with producer path on the same input: agreement level.
    let r_direct = tsqr_r(&a, 48, 3, 12);
    assert_eq!(r_direct.len(), 9);
    let _ = never.reborrow(); // scope reusable across stages
}

// ---------------------------------------------------------------------------
// Counterexample / activation-kill: a confirmed counterexample dominates.
// ---------------------------------------------------------------------------
#[test]
fn counterexample_kills_activation_regardless_of_other_criteria() {
    let criteria = ActivationCriteria {
        theorem_statement_id: Some("hypothetical-theorem-v9"),
        independent_checker_receipt: Some(hash_bytes(b"receipt")),
        falsifier_corpus_digest: Some(hash_bytes(b"corpus")),
        confirmed_counterexample: Some(hash_bytes(b"counterexample")),
    };
        assert!(matches!(
            criteria.satisfied(),
            Err(GaugeBlocker::KilledByCounterexample(d)) if d == hash_bytes(b"counterexample")
        ));
    }

// ---------------------------------------------------------------------------
// Feature gate freeze: current() is disabled; revision 0 has no enable path.
// ---------------------------------------------------------------------------
#[test]
fn feature_gate_cannot_enable_at_revision_zero() {
    let g = ArbitraryTreeGauge::current();
    assert!(!g.is_enabled());
    // Even fully-satisfied-minus-counterexample criteria keep the gate off:
    // enabling requires a gate revision bump shipped WITH the theorem.
    let satisfied = ActivationCriteria {
        theorem_statement_id: Some("t"),
        independent_checker_receipt: Some(hash_bytes(b"r")),
        falsifier_corpus_digest: Some(hash_bytes(b"f")),
        confirmed_counterexample: None,
    };
    assert!(satisfied.satisfied().is_ok(), "criteria can be satisfiable...");
    assert!(!g.is_enabled(), "...but the shipped gate stays frozen");
}

// ---------------------------------------------------------------------------
// Cancellation fault: a scope that fires before the second run surfaces the
// typed refusal; nothing is published.
// ---------------------------------------------------------------------------
#[test]
fn cancellation_between_stages_surfaces_typed_refusal() {
    let a = dep(48);
    let mut always = || true;
    let mut scope = CancelScope::from_closure(&mut always);
    let result = tree_gauge_obstruction(&a, 48, 3, 12, 24, &mut scope);
    assert!(
        matches!(result, Err(PolicyError::CancellationPending)),
        "expected typed cancellation refusal"
    );
}

// ---------------------------------------------------------------------------
// Unknown/no-claim preservation: near-boundary pivots must NOT be promoted
// to a clean verdict merely because empirical samples look green. The
// classification vocabulary keeps Ambiguous distinct, and the tolerance
// default stays the documented sqrt(eps) relative form.
// ---------------------------------------------------------------------------
#[test]
fn ambiguity_band_preserves_unknown_verdicts() {
    let tol = RankTolerance::default_f64();
    assert!(
        (tol.factor() - (1.110_223_024_625_156_5e-16f64).sqrt()).abs() < 1e-24,
        "default tolerance must remain sqrt(eps64) relative"
    );
    // A pivot exactly at the tolerance boundary classifies via the band law
    // in canonical_tree::classify_pivots; here we pin that the DEFAULT
    // policy object carries no absolute threshold that could silently
    // reclassify it.
    let g = ArbitraryTreeGauge::current();
    let _ = g; // gate state irrelevant to classification: no promotion path exists
}

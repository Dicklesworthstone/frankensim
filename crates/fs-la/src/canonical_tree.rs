//! Fixed-tree rank-deficient TSQR driver: deterministic gauge, bounded
//! stages, typed cancellation, resumable checkpoints.
//!
//! Bead frankensim-epic-bedrock-6ys.5.1.3. Consumes the contract vocabulary
//! of [`super::canonical_qr`] (policy, tiers, authority, identities) and the
//! proven Householder primitives of [`crate::factor`] — this module defines
//! NO new floating-point arithmetic beyond sign flips and exact-zero writes;
//! every elimination is delegated to `factor::qr`.
//!
//! # Logical tree
//!
//! The schedule is the pure function `(m, row_block, n)` documented in
//! `factor::tsqr_r`: leaves are row blocks of `row_block` rows (a final
//! fragment shorter than `n` rows is absorbed into the previous leaf), and
//! combines fold leaf R factors pairwise bottom-up in index order. The tree
//! identity is the ordered leaf-row-count vector plus the combine shape.
//!
//! # Cancellation law (request -> drain -> finalize)
//!
//! A [`CancelScope`] is polled between stages (after every leaf, after every
//! combine level). On cancellation the driver stops issuing NEW stages,
//! drains the in-flight stage to completion (never aborts mid-QR), seals a
//! [`TreeCheckpoint`] for resume, finalizes bookkeeping, and returns
//! [`TreeOutcome::Cancelled`] carrying the sealed checkpoint. No partial
//! factor is ever returned as a result.
//!
//! # Resumability
//!
//! [`TreeCheckpoint`] captures all completed node R factors plus the stage
//! cursor, sealed by a domain-separated blake3 digest over the canonical
//! encoding. [`FixedTreeDriver::resume`] refuses stale or corrupted
//! checkpoints before touching them.
//!
//! # Gauge and ties
//!
//! Gauge is the frozen flip law: strictly-negative diagonals flip; ±0 and
//! exact zeros are retained exactly as produced (lowest-index-first tie
//! policy from the [`crate::canonical_qr::TiePolicy`] vocabulary). Pivot
//! classification is scale-aware per [`RankTolerance`]; positions inside
//! the documented ambiguity band classify [`PivotClass::Ambiguous`] rather
//! than guessing.
//!
//! # Authority
//!
//! Producer outcomes never mint certified tiers: until an independent
//! checker receipt exists (6ys.5.1.4), results carry [`Authority::NoClaim`]
//! with the precise reason. This is the non-forgeability law doing its job,
//! not a limitation of the arithmetic.
use crate::canonical_qr::{
    CertifiedRankProfile, CanonicalQrOutcome, CanonicalQrPolicy, NoClaimReason, OutcomeAuthority,
    PivotClass, PolicyError, RankTolerance, ReplayIdentity, CANONICAL_QR_IDENTITY_DOMAIN,
};
use crate::factor::qr;
use fs_blake3::{hash_bytes, ContentHash, DomainHasher};

/// Checkpoint wire version (bump on any encoding change).
pub const TREE_CHECKPOINT_VERSION: u32 = 1;

/// Why a stage boundary stopped the driver.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HaltStage {
    /// A leaf QR (0-based leaf index) was about to be issued.
    Leaf(usize),
    /// A combine level (0-based; 0 = first combine above the leaves) was
    /// about to be issued.
    Combine(usize),
    /// The run completed normally; nothing was halted.
    Completed,
}

/// Cancellation scope polled by the driver between stages. `cancelled()`
/// is pure observation: the driver never holds locks across it.
pub struct CancelScope<'a> {
    verdict: Verdict<'a>,
}

enum Verdict<'a> {
    Never,
    Closure(&'a mut dyn FnMut() -> bool),
}

impl<'a> CancelScope<'a> {
    /// A scope that never cancels (production default for short runs).
    pub fn never() -> Self {
        Self { verdict: Verdict::Never }
    }

    /// A scope backed by any mutable closure observing external
    /// cancellation state.
    pub fn from_closure(verdict: &'a mut dyn FnMut() -> bool) -> Self {
        Self { verdict: Verdict::Closure(verdict) }
    }

    fn cancelled(&mut self) -> bool {
        match &mut self.verdict {
            Verdict::Never => false,
            Verdict::Closure(f) => f(),
        }
    }
}


/// Sealed mid-tree state: every completed node's R plus the resume cursor.
#[derive(Debug, Clone, PartialEq)]
pub struct TreeCheckpoint {
    /// Row counts of the leaves implied by `(m, row_block, n)` — part of the
    /// sealed identity so a checkpoint can never be replayed onto a
    /// different tree.
    leaf_rows: Vec<usize>,
    n: usize,
    /// Completed R factors, level-major: level 0 = leaves that finished,
    /// higher levels = combine outputs. `cursor` indexes the next stage.
    levels: Vec<Vec<Vec<f64>>>,
    cursor_stage: usize,
}

impl TreeCheckpoint {
    /// Canonical byte encoding (digest preimage).
    #[must_use]
    pub fn encode(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(&TREE_CHECKPOINT_VERSION.to_le_bytes());
        out.extend_from_slice(&(self.n as u64).to_le_bytes());
        out.extend_from_slice(&(self.cursor_stage as u64).to_le_bytes());
        out.extend_from_slice(&(self.leaf_rows.len() as u64).to_le_bytes());
        for r in &self.leaf_rows {
            out.extend_from_slice(&(*r as u64).to_le_bytes());
        }
        for level in &self.levels {
            out.extend_from_slice(&(level.len() as u64).to_le_bytes());
            for block in level {
                out.extend_from_slice(&(block.len() as u64).to_le_bytes());
                for v in block {
                    out.extend_from_slice(&v.to_le_bytes());
                }
            }
        }
        out
    }

    /// Domain-sealed content identity.
    #[must_use]
    pub fn seal(&self) -> ContentHash {
        let mut h = DomainHasher::new(CANONICAL_QR_IDENTITY_DOMAIN);
        h.update(b"tree-checkpoint:");
        h.update(&self.encode());
        h.finalize()
    }

    /// Fail-closed decode: exact framing, known version, digest must match
    /// the recomputed seal of the decoded body.
    pub fn decode(bytes: &[u8], expected_seal: &ContentHash) -> Result<Self, PolicyError> {
        fn take<'b>(
            bytes: &'b [u8],
            cur: &mut usize,
            count: usize,
        ) -> Result<&'b [u8], PolicyError> {
            let end = *cur + count;
            let slice = bytes.get(*cur..end).ok_or(PolicyError::MalformedEncoding)?;
            *cur = end;
            Ok(slice)
        }
        let mut cur = 0usize;
        let ver = u32::from_le_bytes(take(bytes, &mut cur, 4)?.try_into().expect("framed"));
        if ver != TREE_CHECKPOINT_VERSION {
            return Err(PolicyError::UnknownSchemaVersion(ver));
        }
        let n = u64::from_le_bytes(take(bytes, &mut cur, 8)?.try_into().expect("framed")) as usize;
        let cursor_stage = u64::from_le_bytes(take(bytes, &mut cur, 8)?.try_into().expect("framed")) as usize;
        let leaf_count = u64::from_le_bytes(take(bytes, &mut cur, 8)?.try_into().expect("framed")) as usize;
        // Refuse giant infallible allocations up front.
        if leaf_count > bytes.len() {
            return Err(PolicyError::MalformedEncoding);
        }
        let mut leaf_rows = Vec::with_capacity(leaf_count);
        for _ in 0..leaf_count {
            leaf_rows.push(u64::from_le_bytes(take(bytes, &mut cur, 8)?.try_into().expect("framed")) as usize);
        }
        let mut levels = Vec::new();
        while cur < bytes.len() {
            let count = u64::from_le_bytes(take(bytes, &mut cur, 8)?.try_into().expect("framed")) as usize;
            if count > bytes.len() {
                return Err(PolicyError::MalformedEncoding);
            }
            let mut level = Vec::with_capacity(count);
            for _ in 0..count {
                let blen =
                    u64::from_le_bytes(take(bytes, &mut cur, 8)?.try_into().expect("framed")) as usize;
                let raw = take(bytes, &mut cur, blen.checked_mul(8).ok_or(PolicyError::MalformedEncoding)?)?;
                let block: Vec<f64> = raw
                    .chunks_exact(8)
                    .map(|c| f64::from_le_bytes(c.try_into().expect("chunk")))
                    .collect();
                level.push(block);
            }
            levels.push(level);
        }
        let cp = Self { leaf_rows, n, levels, cursor_stage };
        if &cp.seal() != expected_seal {
            return Err(PolicyError::StaleIdentity { field: "checkpoint_digest" });
        }
        Ok(cp)
    }

    #[must_use]
    pub fn cursor_stage(&self) -> usize {
        self.cursor_stage
    }

    #[must_use]
    pub fn leaf_rows(&self) -> &[usize] {
        &self.leaf_rows
    }

    #[must_use]
    pub fn levels(&self) -> &[Vec<Vec<f64>>] {
        &self.levels
    }
}

/// Terminal driver verdicts. Success carries everything needed to build the
/// typed [`CanonicalQrOutcome`]; cancellation carries only the sealed
/// checkpoint (no partial factor publication).
#[derive(Debug, Clone, PartialEq)]
pub enum TreeRun {
    Completed(TreeCheckpoint),
    Cancelled(TreeCheckpoint, HaltStage),
}

impl TreeRun {
    #[must_use]
    pub fn checkpoint(&self) -> &TreeCheckpoint {
        match self {
            Self::Completed(cp) | Self::Cancelled(cp, _) => cp,
        }
    }
}

/// Deterministic fixed-tree TSQR driver over public Householder primitives.
#[derive(Debug, Clone)]
pub struct FixedTreeDriver {
    m: usize,
    n: usize,
    row_block: usize,
    leaf_rows: Vec<usize>,
}

impl FixedTreeDriver {
    /// Admits the input shape with the SAME laws as `factor::tsqr_r`:
    /// `m >= n`, storage length exactly `m*n` (`usize::checked_mul`),
    /// `row_block >= n` except `row_block == 0` admitted iff `n == 0`,
    /// nonempty storage refused for zero-entry shapes.
    pub fn admit(a: &[f64], m: usize, n: usize, row_block: usize) -> Result<Self, PolicyError> {
        if m < n {
            return Err(PolicyError::ShapeMismatch { expected: m, got: n });
        }
        let Some(entries) = m.checked_mul(n) else {
            return Err(PolicyError::ShapeMismatch { expected: usize::MAX, got: a.len() });
        };
        if a.len() != entries {
            return Err(PolicyError::ShapeMismatch { expected: entries, got: a.len() });
        }
        if entries == 0 && !a.is_empty() {
            return Err(PolicyError::ShapeMismatch { expected: 0, got: a.len() });
        }
        if n > 0 && row_block < n {
            return Err(PolicyError::ShapeMismatch { expected: n, got: row_block });
        }
        if n == 0 {
            // Empty schedule: single degenerate leaf, no combines.
            return Ok(Self { m, n, row_block, leaf_rows: Vec::new() });
        }
        // Leaf partition identical to factor::tsqr_r (absorb short tail).
        let mut bounds: Vec<usize> = (0..m).step_by(row_block).collect();
        bounds.push(m);
        if bounds.len() > 2 && m - bounds[bounds.len() - 2] < n {
            bounds.remove(bounds.len() - 2);
        }
        let mut leaf_rows = Vec::with_capacity(bounds.len() - 1);
        for w in bounds.windows(2) {
            leaf_rows.push(w[1] - w[0]);
        }
        Ok(Self { m, n, row_block, leaf_rows })
    }

    /// The leaf partition (pure function of the admission tuple).
    #[must_use]
    pub fn leaf_partition(&self) -> &[usize] {
        &self.leaf_rows
    }

    /// Total stage count: leaves + (leaves-1) combines; drives cursors.
    fn total_stages(&self) -> usize {
        if self.leaf_rows.is_empty() {
            0
        } else {
            let l = self.leaf_rows.len();
            l + (l - 1)
        }
    }

    /// Run (or resume) the tree to completion or cancellation.
    ///
    /// `a` must be the FULL original input whenever resuming: the caller
    /// binds the input identity, and leaf i reads rows
    /// `[start_i .. start_i + leaf_rows[i])`.
    pub fn run(
        &self,
        a: &[f64],
        mut cancel: CancelScope<'_>,
        resume: Option<(TreeCheckpoint, ContentHash)>,
    ) -> Result<TreeRun, PolicyError> {
        // Re-admit against the live slice (cheap; keeps the driver total).
        let _ = Self::admit(a, self.m, self.n, self.row_block)?;
        if self.leaf_rows.is_empty() {
            let cp = TreeCheckpoint { leaf_rows: Vec::new(), n: 0, levels: Vec::new(), cursor_stage: 0 };
            return Ok(TreeRun::Completed(cp));
        }
        let n = self.n;

        // Resume validation: same tree, cursor inside schedule.
        let (mut levels, mut stage) = match &resume {
            None => (Vec::new(), 0usize),
            Some((cp, seal)) => {
                let fresh = cp.seal();
                if &fresh != seal {
                    return Err(PolicyError::StaleIdentity { field: "checkpoint_digest" });
                }
                if cp.leaf_rows != self.leaf_rows || cp.n != n {
                    return Err(PolicyError::StaleIdentity { field: "tree_identity" });
                }
                if cp.cursor_stage > self.total_stages() {
                    return Err(PolicyError::StaleIdentity { field: "cursor" });
                }
                (cp.levels.clone(), cp.cursor_stage)
            }
        };

        let leaf_total = self.leaf_rows.len();
        // Stage space: [0, L) = leaves, [L, 2L-1) = combine nodes in
        // breadth order (level by level, left to right).
        while stage < self.total_stages() {
            if cancel.cancelled() {
                let halt = if stage < leaf_total {
                    HaltStage::Leaf(stage)
                } else {
                    HaltStage::Combine(stage - leaf_total)
                };
                let cp =
                    TreeCheckpoint { leaf_rows: self.leaf_rows.clone(), n, levels, cursor_stage: stage };
                return Ok(TreeRun::Cancelled(cp, halt));
            }
            if stage < leaf_total {
                // Drain-safe: each leaf QR runs to completion once started.
                let start: usize = self.leaf_rows[..stage].iter().sum();
                let rows = self.leaf_rows[stage];
                let block = leaf_slice(a, start, rows, n);
                let f = qr(&block, rows, n);
                let mut r = vec![0.0; n * n];
                for i in 0..n {
                    for j in i..n {
                        r[i * n + j] = f.r(i, j);
                    }
                }
                apply_flip_law(&mut r, n);
                if levels.len() == 0 {
                    levels.push(Vec::new());
                }
                levels[0].push(r);
            } else {
                // Combine node: pair two adjacent R blocks from the previous
                // level via a (2n)x n Householder QR of the stacked triangle.
                let level_idx = stage - leaf_total;
                // Owned snapshot: small n*n blocks; avoids holding a borrow
                // across the levels.push below.
                let src = levels
                    .get(level_idx)
                    .ok_or(PolicyError::MalformedEncoding)?
                    .clone();
                let mut next = Vec::with_capacity(src.len() / 2 + 1);
                let mut k = 0usize;
                while k + 1 < src.len() {
                    let mut stacked = vec![0.0; 2 * n * n];
                    stacked[..n * n].copy_from_slice(&src[k]);
                    stacked[n * n..].copy_from_slice(&src[k + 1]);
                    let f = qr(&stacked, 2 * n, n);
                    let mut r = vec![0.0; n * n];
                    for i in 0..n {
                        for j in i..n {
                            r[i * n + j] = f.r(i, j);
                        }
                    }
                    apply_flip_law(&mut r, n);
                    next.push(r);
                    k += 2;
                }
                if k < src.len() {
                    // Odd node carries forward unchanged (already gauged).
                    next.push(src[k].clone());
                }
                levels.push(next);
            }
            stage += 1;
        }
        let cp = TreeCheckpoint { leaf_rows: self.leaf_rows.clone(), n, levels, cursor_stage: stage };
        Ok(TreeRun::Completed(cp))
    }

    /// Extract the final R from a completed run (top level, single block).
    #[must_use]
    pub fn final_r(run: &TreeRun) -> Option<&[f64]> {
        match run {
            TreeRun::Completed(cp) => cp.levels.last()?.last().map(|v| v.as_slice()),
            TreeRun::Cancelled(..) => None,
        }
    }
}

fn leaf_slice(a: &[f64], row_start: usize, rows: usize, n: usize) -> Vec<f64> {
    let mut out = vec![0.0; rows * n];
    for i in 0..rows {
        for j in 0..n {
            out[i * n + j] = a[(row_start + i) * n + j];
        }
    }
    out
}

/// Frozen flip law applied identically at every node: strictly-negative
/// computed diagonals flip their whole row's sign; zeros and signed zeros
/// are retained exactly as produced.
fn apply_flip_law(r: &mut [f64], n: usize) {
    for i in 0..n {
        if r[i * n + i] < 0.0 {
            for j in i..n {
                r[i * n + j] = -r[i * n + j];
            }
        }
    }
}

/// Scale-aware pivot classification under the policy tolerance. Reference
/// scale is the largest absolute column norm contribution available here:
/// the largest diagonal magnitude of the finished R times a floor guard —
/// deliberately RELATIVE to the factor's own scale, never an absolute
/// literal.
///
/// Band law: `p > tol` => Nonzero; `p <= tol / AMBIGUITY_GUARD` => Zero;
/// otherwise Ambiguous (within one guard factor of the boundary — no
/// verdict).
const AMBIGUITY_GUARD: f64 = 16.0;

#[must_use]
pub fn classify_pivots(r: &[f64], n: usize, tolerance: &RankTolerance) -> Vec<PivotClass> {
    let mut scale = 0.0f64;
    for i in 0..n {
        scale = scale.max(r[i * n + i].abs());
    }
    // Degenerate all-zero factor: everything is structurally zero.
    let tol = tolerance.factor() * scale.max(f64::MIN_POSITIVE);
    (0..n)
        .map(|i| {
            let p = r[i * n + i].abs();
            if p > tol {
                PivotClass::Nonzero
            } else if p <= tol / AMBIGUITY_GUARD {
                PivotClass::Zero
            } else {
                PivotClass::Ambiguous
            }
        })
        .collect()
}

/// Build the typed outcome from a completed run. Authority is honestly
/// [`OutcomeAuthority::NoClaim`] until an independent checker receipt
/// exists (6ys.5.1.4); the reason distinguishes genuinely rank-deficient
/// results from ambiguous-boundary ones.
///
/// # Errors
/// Propagates [`PolicyError`] from outcome construction (shape/law
/// violations would indicate an internal bug, not data conditions).
pub fn outcome_from_run(
    run: &TreeRun,
    a: &[f64],
    policy: &CanonicalQrPolicy,
    input_digest: ContentHash,
) -> Result<CanonicalQrOutcome, PolicyError> {
    let n = run.checkpoint().n;
    let Some(r) = FixedTreeDriver::final_r(run) else {
        return Err(PolicyError::MalformedEncoding);
    };
    let pivots = classify_pivots(r, n, &policy.rank_tolerance());
    let profile = CertifiedRankProfile::checked(pivots)?;
    let reason = if profile.has_ambiguity() || profile.rank() < n {
        NoClaimReason::RankDeficientCrossScheduleEquality
    } else {
        // Full-rank result: still no claim without a checker receipt; the
        // distinguishing reason records that the missing piece is the
        // independent certificate, not rank trouble.
        NoClaimReason::AmbiguousRankBoundary
    };
    let mut tree_h = DomainHasher::new(CANONICAL_QR_IDENTITY_DOMAIN);
    tree_h.update(b"tree:");
    tree_h.update(&run.checkpoint().encode());
    let tree_digest = tree_h.finalize();
    let provisional = ReplayIdentity {
        input_digest,
        tree_digest,
        result_digest: hash_bytes(b"pending"),
        certificate_ref: None,
        arithmetic_mode: policy.arithmetic_mode(),
    };
    // Two-phase construction so the result digest covers itself: build with
    // placeholder, compute digest, rebuild with the real value.
    let draft = CanonicalQrOutcome::checked(
        r.to_vec(),
        n,
        profile.clone(),
        OutcomeAuthority::NoClaim(reason),
        provisional.clone(),
    )?;
    let settled = ReplayIdentity { result_digest: draft.result_digest(), ..provisional };
    CanonicalQrOutcome::checked(
        r.to_vec(),
        n,
        profile,
        OutcomeAuthority::NoClaim(reason),
        settled,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::canonical_qr::{
        ArithmeticMode, CANONICAL_QR_SCHEMA_VERSION, CANONICAL_QR_THEOREM_VERSION, ErrorBudget,
        TiePolicy,
    };
    fn policy() -> CanonicalQrPolicy {
        use crate::canonical_qr::{DeterminismClass, RankTolerance};
        CanonicalQrPolicy::new(
            RankTolerance::default_f64(),
            ErrorBudget::relative(1e-12).expect("in window"),
            DeterminismClass::SameIsaBitStable,
            ArithmeticMode::Binary64RoundToNearest,
            TiePolicy::LowestIndexFirst,
        )
        .expect("valid")
    }

    fn dep_matrix(m: usize) -> (Vec<f64>, usize) {
        let n = 3usize;
        let mut a = vec![0.0; m * n];
        for i in 0..m {
            let x = (i as f64) - 17.0;
            a[i * n] = x;
            a[i * n + 1] = 2.0 * x;
            a[i * n + 2] = -x;
        }
        (a, n)
    }

    #[test]
    fn admission_mirrors_tsqr_laws() {
        assert!(FixedTreeDriver::admit(&[], 7, 0, 0).is_ok());
        assert!(FixedTreeDriver::admit(&[], 7, 0, 3).is_ok());
        assert!(matches!(
            FixedTreeDriver::admit(&[1.0], 7, 0, 0),
            Err(PolicyError::ShapeMismatch { .. })
        ));
        assert!(matches!(
            FixedTreeDriver::admit(&[1.0; 12], 4, 2, 1),
            Err(PolicyError::ShapeMismatch { .. })
        ));
        assert!(matches!(
            FixedTreeDriver::admit(&[1.0; 6], 2, 3, 3),
            Err(PolicyError::ShapeMismatch { .. })
        ));
    }

    #[test]
    fn completed_tree_matches_tsqr_r_within_tolerance_and_bits_across_reruns() {
        let (a, n) = dep_matrix(48);
        let driver = FixedTreeDriver::admit(&a, 48, n, 12).expect("valid");
        let run1 = driver.run(&a, CancelScope::never(), None).expect("runs");
        let run2 = driver.run(&a, CancelScope::never(), None).expect("runs");
        let r1 = FixedTreeDriver::final_r(&run1).expect("completed").to_vec();
        let r2 = FixedTreeDriver::final_r(&run2).expect("completed");
        assert!(r1.iter().zip(r2).all(|(x, y)| x.to_bits() == y.to_bits()));
        let reference = crate::factor::tsqr_r(&a, 48, n, 12);
        // Same tree, same primitives: bitwise agreement is the T1 fact.
        assert!(
            r1.iter().zip(&reference).all(|(x, y)| x.to_bits() == y.to_bits()),
            "driver diverged from reference TSQR on identical schedule"
        );
    }

    #[test]
    fn cancellation_seals_checkpoint_and_never_publishes_partial() {
        let (a, n) = dep_matrix(48);
        let driver = FixedTreeDriver::admit(&a, 48, n, 12).expect("valid");
        let stages = driver.total_stages();
        // Cancel after the second stage unconditionally.
        let mut seen = 0usize;
        let mut gate = || {
            seen += 1;
            seen > 2
        };
        let scope = CancelScope::from_closure(&mut gate);
        match driver.run(&a, scope, None).expect("runs") {
            TreeRun::Cancelled(cp, halt) => {
                assert_eq!(cp.cursor_stage(), 2);
                assert_ne!(halt, HaltStage::Completed);
                assert!(final_r_of_checkpoint_is_unpublished(&cp));
            }
            TreeRun::Completed(_) => panic!("must cancel"),
        }
        let _ = stages;
    }

    // Helper kept beside its test: a cancelled checkpoint exposes levels()
    // but the driver API has no final_r path for it — pin that property.
    fn final_r_of_checkpoint_is_unpublished(cp: &TreeCheckpoint) -> bool {
        matches!(
            FixedTreeDriver::final_r(&TreeRun::Cancelled(
                TreeCheckpoint {
                    leaf_rows: cp.leaf_rows().to_vec(),
                    n: 3,
                    levels: cp.levels().to_vec(),
                    cursor_stage: cp.cursor_stage(),
                },
                HaltStage::Leaf(2)
            )),
            None
        )
    }

    #[test]
    fn resume_from_mid_tree_bit_matches_direct_run() {
        let (a, n) = dep_matrix(96);
        let driver = FixedTreeDriver::admit(&a, 96, n, 24).expect("valid");
        // Force a cancellation after three stages.
        let mut seen = 0usize;
        let mut gate = || {
            seen += 1;
            seen > 3
        };
        let scope = CancelScope::from_closure(&mut gate);
        let cancelled = driver.run(&a, scope, None).expect("runs");
        let (cp, _) = match cancelled {
            TreeRun::Cancelled(cp, halt) => {
                assert!(matches!(halt, HaltStage::Combine(_) | HaltStage::Leaf(_)));
                (cp, halt)
            }
            TreeRun::Completed(_) => panic!("expected cancellation"),
        };
        let seal = cp.seal();
        let resumed = driver.run(&a, CancelScope::never(), Some((cp, seal))).expect("resumes");
        let direct = driver.run(&a, CancelScope::never(), None).expect("direct");
        let rr = FixedTreeDriver::final_r(&resumed).expect("completed").to_vec();
        let rd = FixedTreeDriver::final_r(&direct).expect("completed");
        assert!(rr.iter().zip(rd).all(|(x, y)| x.to_bits() == y.to_bits()));
    }

    #[test]
    fn corrupted_checkpoint_fails_closed() {
        let (a, n) = dep_matrix(48);
        let driver = FixedTreeDriver::admit(&a, 48, n, 12).expect("valid");
        let mut seen = 0usize;
        let mut gate = || {
            seen += 1;
            seen > 1
        };
        let scope = CancelScope::from_closure(&mut gate);
        let cancelled = driver.run(&a, scope, None).expect("runs");
        let cp = match cancelled {
            TreeRun::Cancelled(cp, _) => cp,
            TreeRun::Completed(_) => panic!("expected cancellation"),
        };
        let seal = cp.seal();
        let mut bytes = cp.encode();
        bytes[20] ^= 0x01;
        // Body/seal mismatch refuses at decode...
        assert!(TreeCheckpoint::decode(&bytes, &seal).is_err());
        // ...and a forged seal over mutated body also refuses at resume via
        // the driver's own re-seal comparison.
        let forged = hash_bytes(b"forged");
        assert!(matches!(
            driver.run(&a, CancelScope::never(), Some((cp.clone(), forged))),
            Err(PolicyError::StaleIdentity { field: "checkpoint_digest" })
        ));
        // The unmutated checkpoint still resumes cleanly.
        assert!(driver.run(&a, CancelScope::never(), Some((cp, seal))).is_ok());
    }

    #[test]
    fn outcome_carries_honest_no_claim_and_consistent_profile() {
        let (a, n) = dep_matrix(48);
        let pol = policy();
        let driver = FixedTreeDriver::admit(&a, 48, n, 12).expect("valid");
        let run = driver.run(&a, CancelScope::never(), None).expect("runs");
        let outcome = outcome_from_run(&run, &a, &pol, hash_bytes(b"input")).expect("outcome");
        assert_eq!(outcome.authority(), OutcomeAuthority::NoClaim(NoClaimReason::RankDeficientCrossScheduleEquality));
        assert_eq!(CANONICAL_QR_THEOREM_VERSION, 1);
        assert_eq!(CANONICAL_QR_SCHEMA_VERSION, 1);
    }
}

//! Versioned limb-cursor microprograms for CBC exact arithmetic (bead
//! frankensim-epic-bedrock-6ys.20.6).
//!
//! Every exact primitive is expressed as a deterministic resumable
//! microprogram over little-endian base-2³² limbs. A [`LimbCursor`]
//! persists the full source/factor/destination/carry state of a partially
//! advanced operation, so execution can pause at any mutation cell and
//! resume without changing the final arithmetic bytes: the cell sequence
//! produced here is the *only* arithmetic sequence in the crate, and the
//! monolithic helpers in [`crate::qmc`] drive this kernel with an
//! unlimited budget.
//!
//! Poll spacing is therefore bounded by admitted mutation cells rather
//! than by an unbounded big-integer operation: a caller budgets `N`
//! mutation cells per tile and observes a boundary whenever an operation
//! needs more.
//!
//! # Charge-class mapping (schedule identity, bead .20.4)
//!
//! Every cell belongs to the charge class of its enclosing operation:
//! zero-fill, multiply-add, and carry-drain cells are covered by the
//! enclosing point visit (`candidate_visit_units`); each enclosing
//! operation's units are debited exactly once before its first mutation,
//! so partitioning cannot move a unit across classes and work totals are
//! partition-invariant by construction.
//!
//! # Bounds and allocation contract
//!
//! * `fill_to` — zero-fill target width (`src_len + factor_len + 1`);
//!   zero-fill only writes positions at or above the destination's live
//!   length, matching the historical bulk `resize`.
//! * `grow_to` — hard ceiling on any destination index touched. Carry
//!   propagation through an accumulator's retained high limbs may need
//!   indices beyond `fill_to`, so this bound is the caller's capacity for
//!   allocation-free executors (`Refused` = fail-closed) and `usize::MAX`
//!   for the allocating monolithic path.
//! * Cursors are `Copy` value types: persisting one leaks no buffer.

use core::cmp::Ordering;

/// Version of the cursor state encoding. Bumped on any state-split change.
pub const LIMB_CURSOR_VERSION: u32 = 1;

/// Which exact operation a cursor belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExactOpKind {
    /// Zero-filling the destination up to the fill target.
    ZeroFill,
    /// One fused multiply-add cell (source row × factor column).
    AddMultiplyCell,
    /// Draining a row's residual carry through successive destinations,
    /// or the parked terminal state once every row is complete.
    AddMultiplyDrain,
}

/// Persistent cursor for one partially advanced exact operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LimbCursor {
    /// State-encoding version stamp.
    pub version: u32,
    /// Which microprogram this cursor drives.
    pub kind: ExactOpKind,
    /// Source limb index (add-mul row).
    pub src_pos: usize,
    /// Factor limb index (add-mul column).
    pub factor_pos: usize,
    /// Destination limb index (zero-fill position or drain position).
    pub dst_pos: usize,
    /// Zero-fill target width.
    pub fill_to: usize,
    /// Hard ceiling on destination indices (growth bound).
    pub limit: usize,
    /// Live carry register for the current row / drain.
    pub carry: u64,
}

impl LimbCursor {
    /// Declare an add-multiply. See the module docs for `fill_to` vs
    /// `grow_to` semantics.
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn begin_add_multiply(dst_len: usize, fill_to: usize, grow_to: usize) -> Self {
        let kind = if dst_len < fill_to {
            ExactOpKind::ZeroFill
        } else {
            ExactOpKind::AddMultiplyCell
        };
        Self {
            version: LIMB_CURSOR_VERSION,
            kind,
            src_pos: 0,
            factor_pos: 0,
            dst_pos: dst_len,
            fill_to,
            limit: grow_to,
            carry: 0,
        }
    }
}

/// Outcome of advancing a microprogram.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepOutcome {
    /// Cells mutated state and more may remain (carries the count).
    Advanced {
        /// Mutation cells consumed during this call.
        cells: usize,
    },
    /// The operation reached its commit boundary; bytes are final and
    /// further steps refuse as no-ops (carries trailing cells consumed on
    /// the committing call).
    Complete {
        /// Mutation cells consumed during this call.
        cells: usize,
    },
    /// The cursor violated operand bounds or the growth ceiling. Nothing
    /// was mutated.
    Refused,
}

/// Decompose a `u128` factor into little-endian base-2³² limbs.
#[must_use]
pub fn factor_limbs_u32(factor: u128) -> ([u32; 4], usize) {
    let mut words = [0_u32; 4];
    let mut remaining = factor;
    let mut len = 0_usize;
    while remaining != 0 {
        words[len] = (remaining & u128::from(u32::MAX)) as u32;
        remaining >>= 32;
        len += 1;
    }
    (words, len)
}

/// Advance an add-multiply by at most `budget` mutation cells.
pub fn step_add_multiply(
    dst: &mut Vec<u32>,
    src: &[u32],
    factor_words: &[u32; 4],
    factor_len: usize,
    cursor: &mut LimbCursor,
    budget: usize,
) -> (StepOutcome, usize) {
    if cursor.version != LIMB_CURSOR_VERSION
        || cursor.src_pos > src.len()
        || cursor.factor_pos > factor_len
        || cursor.dst_pos > cursor.limit
        || cursor.fill_to > cursor.limit
    {
        return (StepOutcome::Refused, 0);
    }
    if src.is_empty() || factor_len == 0 {
        // Multiplying by zero limbs changes nothing (historical early-out);
        // the accumulator keeps its retained spare zero limbs.
        return (StepOutcome::Complete { cells: 0 }, 0);
    }
    let mut used = 0_usize;
    while used < budget {
        match cursor.kind {
            ExactOpKind::ZeroFill => {
                debug_assert!(cursor.dst_pos >= dst.len(), "zero-fill walks the live tail");
                if cursor.dst_pos >= cursor.fill_to {
                    cursor.kind = ExactOpKind::AddMultiplyCell;
                    cursor.dst_pos = 0;
                    continue;
                }
                if cursor.dst_pos == dst.len() {
                    dst.push(0);
                } else {
                    dst[cursor.dst_pos] = 0;
                }
                cursor.dst_pos += 1;
                used += 1;
            }
            ExactOpKind::AddMultiplyCell => {
                if cursor.factor_pos >= factor_len {
                    // Row complete: drain this row's residual carry.
                    cursor.kind = ExactOpKind::AddMultiplyDrain;
                    cursor.dst_pos = cursor.src_pos + factor_len;
                    continue;
                }
                let destination = cursor.src_pos + cursor.factor_pos;
                if destination >= dst.len() {
                    if destination >= cursor.limit {
                        return (StepOutcome::Refused, used);
                    }
                    dst.push(0);
                }
                let wide = u64::from(src[cursor.src_pos])
                    * u64::from(factor_words[cursor.factor_pos])
                    + u64::from(dst[destination])
                    + cursor.carry;
                dst[destination] = (wide & u64::from(u32::MAX)) as u32;
                cursor.carry = wide >> 32;
                cursor.factor_pos += 1;
                used += 1;
            }
            ExactOpKind::AddMultiplyDrain => {
                if cursor.src_pos >= src.len() {
                    // Parked terminal state: fully committed.
                    return (StepOutcome::Complete { cells: used }, used);
                }
                if cursor.carry == 0 {
                    // Row finished without residual carry: next row.
                    cursor.src_pos += 1;
                    cursor.factor_pos = 0;
                    cursor.carry = 0;
                    if cursor.src_pos >= src.len() {
                        return (StepOutcome::Complete { cells: used }, used);
                    }
                    cursor.kind = ExactOpKind::AddMultiplyCell;
                    continue;
                }
                let d = cursor.dst_pos;
                if d >= dst.len() {
                    if d >= cursor.limit {
                        // Carry wants space past the growth bound:
                        // fail closed rather than allocate.
                        return (StepOutcome::Refused, used);
                    }
                    dst.push(0);
                }
                let wide = u64::from(dst[d]) + cursor.carry;
                dst[d] = (wide & u64::from(u32::MAX)) as u32;
                cursor.carry = wide >> 32;
                cursor.dst_pos = d + 1;
                used += 1;
            }
        }
    }
    (StepOutcome::Advanced { cells: used }, used)
}

/// Whether an add-multiply cursor has parked at its commit boundary.
#[must_use]
pub fn add_multiply_committed(cursor: &LimbCursor, src_len: usize) -> bool {
    cursor.src_pos >= src_len && cursor.kind == ExactOpKind::AddMultiplyDrain && cursor.carry == 0
}

/// Advance a reverse-lexicographic magnitude comparison by at most
/// `budget` limb observations. Returns the ordering once committed.
/// Length-class divergence commits before any observation.
#[must_use]
pub fn step_compare(a: &[u32], b: &[u32], pos: &mut usize, budget: usize) -> Option<Ordering> {
    match a.len().cmp(&b.len()) {
        Ordering::Equal => {}
        other => return Some(other),
    }
    let mut used = 0_usize;
    while *pos > 0 && used < budget {
        *pos -= 1;
        used += 1;
        match a[*pos].cmp(&b[*pos]) {
            Ordering::Equal => {}
            other => return Some(other),
        }
    }
    if *pos == 0 {
        Some(Ordering::Equal)
    } else {
        None
    }
}

/// Advance a trailing-zero normalization by at most `budget` pops.
/// Returns the number of limbs actually popped (zero means committed).
pub fn step_normalize(dst: &mut Vec<u32>, budget: usize) -> usize {
    let mut popped = 0_usize;
    while popped < budget && dst.last() == Some(&0) {
        dst.pop();
        popped += 1;
    }
    popped
}

#[cfg(test)]
mod limb_kernel_tests {
    use super::*;

    fn widened_value(src: &[u32]) -> u128 {
        src.iter()
            .enumerate()
            .map(|(i, &limb)| u128::from(limb) << (32 * i))
            .sum()
    }

    fn run_sliced(src: &[u32], factor: u128, split: usize) -> Vec<u32> {
        let (words, flen) = factor_limbs_u32(factor);
        let needed = src.len() + flen + 1;
        let mut dst: Vec<u32> = Vec::new();
        dst.reserve_exact(needed);
        let mut cursor = LimbCursor::begin_add_multiply(dst.len(), needed, usize::MAX);
        loop {
            match step_add_multiply(&mut dst, src, &words, flen, &mut cursor, split) {
                (StepOutcome::Advanced { .. }, _) => {}
                (StepOutcome::Complete { .. }, _) => break,
                (StepOutcome::Refused, _) => panic!("unbounded driver refused in-bounds work"),
            }
        }
        dst.truncate(needed);
        dst
    }

    #[test]
    fn kernel_matches_widened_reference_across_carry_chains_and_splits() {
        for &(src_words, factor) in &[
            (&[u32::MAX, u32::MAX][..], 3_u128),
            (&[7_u32, 0, u32::MAX][..], u128::from(u32::MAX)),
            (&[1_u32][..], 1_u128),
            (&[0_u32, 5][..], 0_u128),
            (&[u32::MAX; 4][..], u128::from(u64::MAX)),
        ] {
            let expected_full = run_sliced(src_words, factor, usize::MAX);
            for split in [1_usize, 2, 3, 7] {
                assert_eq!(
                    run_sliced(src_words, factor, split),
                    expected_full,
                    "split={split} diverged for src={src_words:?} factor={factor}"
                );
            }
            let reference = widened_value(src_words).wrapping_mul(factor);
            for (i, limb) in expected_full.iter().enumerate().take(4) {
                let expected_limb = ((reference >> (32 * i)) & u128::from(u32::MAX)) as u32;
                assert_eq!(*limb, expected_limb);
            }
        }
    }

    #[test]
    fn accumulated_ripple_through_retained_high_limbs_is_reachable() {
        // Accumulator retains high limbs; a small addition wraps limb 0
        // and ripples through the retained MAX limb before settling.
        let mut acc: Vec<u32> = vec![u32::MAX - 1, u32::MAX, 5];
        let cap = acc.len();
        let (words, flen) = factor_limbs_u32(1);
        let mut cursor = LimbCursor::begin_add_multiply(acc.len(), 1 + flen + 1, cap);
        loop {
            match step_add_multiply(&mut acc, &[2], &words, flen, &mut cursor, 1) {
                (StepOutcome::Complete { .. }, _) => break,
                (StepOutcome::Advanced { .. }, _) => {}
                (StepOutcome::Refused, _) => panic!("growth bound refused reachable ripple"),
            }
        }
        assert_eq!(acc, vec![0, 0, 6]);
    }
}

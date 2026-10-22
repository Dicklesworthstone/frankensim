//! Portable checkpoints for feasibility-first differential evolution.
//!
//! The shared DE codec retains exact candidates, finite scores, spent attempts,
//! parent/replacement populations and successful siblings of interrupted waves.
//! This format additionally binds the row count and exact tolerance bits. The
//! feasible flag is reconstructed from the retained violation and that policy;
//! it cannot be imported independently or changed to bypass a constraint.
//!
//! Model identity must cover the objective, all inequality definitions/scales,
//! physical inputs, and numerical implementation. The caller owns that identity
//! and crash-safe storage. A checksum detects corruption, not malicious edits;
//! snapshots do not attest that stored scores were produced by actual physics.

use super::{ConstrainedDifferentialEvolutionState, FeasibilityScore,
    DifferentialEvolutionScore, RankedDifferentialEvolutionState, Error, invalid};
use super::super::snapshot::{SnapshotCodec, Encoder, Decoder};
use fs_exec::Cx;

pub use super::super::snapshot::DifferentialEvolutionRestoreLimits;

/// Independent format version; scalar DE snapshot bytes are unchanged.
pub const CONSTRAINED_DIFFERENTIAL_EVOLUTION_CHECKPOINT_VERSION: u32 = 1;
const MAGIC: &[u8; 8] = b"FSDECF01";
const DOMAIN: &str = "fs-dfo/constrained-differential-evolution-checkpoint-v1";

impl ConstrainedDifferentialEvolutionState {
    /// Immutable declarations for binding a reconstructed physical owner.
    pub fn configuration(&self) -> &super::DifferentialEvolutionConfig { &self.state.config }
    pub fn bounds(&self) -> &[[f64; 2]] { &self.state.bounds }
    pub fn initial_coordinates(&self) -> Option<&[f64]> { self.state.initial.as_deref() }

    /// Capture only a drained search. No callback or budget charge occurs.
    /// Successful cached siblings, failed-candidate retry coordinates and all
    /// charged attempts remain available after restoring in a fresh process.
    /// After cancellation, use a fresh noncancelled Cx for capture.
    pub fn checkpoint_bytes(&self, cx: &Cx<'_>, model_identity: [u8; 32], max_bytes: usize)
        -> Result<Vec<u8>, Error>
    {
        let codec = FeasibilityCodec::new(self.constraint_count, self.feasibility_tolerance)?;
        self.state.encode_snapshot_with(model_identity, max_bytes, &codec,
            &mut || cx.checkpoint().is_ok())
    }

    /// Restore without evaluating physics or replenishing spent work. Both
    /// model identity and caller-supplied feasibility policy must match exactly.
    /// Byte/shape limits are supplied independently of the checkpoint. Unknown
    /// versions, malformed scores/counters, changed policies and foreign cached
    /// candidates refuse before a state is returned. Scalar snapshots refuse.
    ///
    /// A saved exhausted allowance stays exhausted. Explicitly extend it before
    /// requesting more callbacks. Restoring an older snapshot cannot account
    /// for work after capture; a cross-process work ledger is application-owned.
    pub fn from_checkpoint_bytes(cx: &Cx<'_>, bytes: &[u8], expected_model: [u8; 32],
        expected_constraint_count: usize, expected_feasibility_tolerance: f64,
        limits: DifferentialEvolutionRestoreLimits) -> Result<Self, Error>
    {
        let codec = FeasibilityCodec::new(expected_constraint_count, expected_feasibility_tolerance)?;
        let state = RankedDifferentialEvolutionState::<FeasibilityScore>::decode_snapshot_with(
            bytes, expected_model, limits, &codec, &mut || cx.checkpoint().is_ok())?;
        Ok(Self { state, constraint_count: codec.count, feasibility_tolerance: codec.tolerance })
    }
}

struct FeasibilityCodec { count: usize, tolerance: f64 }
impl FeasibilityCodec {
    fn new(count: usize, tolerance: f64) -> Result<Self, Error> {
        if count == 0 || !tolerance.is_finite() || tolerance < 0.0 {
            return Err(invalid("constrained DE checkpoint requires a positive row count and finite nonnegative tolerance"));
        }
        Ok(Self { count, tolerance })
    }
}
impl SnapshotCodec<FeasibilityScore> for FeasibilityCodec {
    const MAGIC: &'static [u8; 8] = self::MAGIC;
    const VERSION: u32 = CONSTRAINED_DIFFERENTIAL_EVOLUTION_CHECKPOINT_VERSION;
    const DOMAIN: &'static str = self::DOMAIN;

    fn write_policy(&self, out: &mut Encoder) -> Result<(), Error> {
        out.count(self.count)?;
        out.number(self.tolerance)
    }
    fn read_policy(&self, input: &mut Decoder<'_>) -> Result<(), Error> {
        if input.count()? != self.count || input.number()?.to_bits() != self.tolerance.to_bits() {
            return Err(invalid("constrained DE checkpoint row count or feasibility tolerance changed"));
        }
        Ok(())
    }
    fn write_score(&self, out: &mut Encoder, score: FeasibilityScore) -> Result<(), Error> {
        let score = score.checked()?;
        if score.feasible != (score.maximum_violation <= self.tolerance) {
            return Err(invalid("constrained DE checkpoint score disagrees with its feasibility policy"));
        }
        out.number(score.objective)?;
        out.number(score.maximum_violation)
    }
    fn read_score(&self, input: &mut Decoder<'_>) -> Result<FeasibilityScore, Error> {
        let objective = input.number()?;
        let maximum_violation = input.number()?;
        FeasibilityScore { objective, maximum_violation,
            feasible: maximum_violation <= self.tolerance }.checked()
    }
    fn equivalent_best(&self, best: FeasibilityScore, minimum: FeasibilityScore) -> bool {
        // The earliest feasible best may no longer be a population member:
        // equal-objective feasible trials replace parents even if their raw
        // violations differ. Validate rank, not objective-minimum or all bits.
        best.rank(minimum).is_eq()
    }
}

#[cfg(test)]
mod tests;

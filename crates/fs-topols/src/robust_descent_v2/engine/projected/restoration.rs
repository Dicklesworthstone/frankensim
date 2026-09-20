//! Explicit sampled-stress feasibility restoration before compliance descent.
//!
//! While infeasible, the governing stress case supplies an unweighted
//! compliance-based proposal. Acceptance uses the complete re-solved family's
//! ACTUAL sampled maximum, not a derivative prediction. This heuristic can
//! stall; it is not a stress adjoint or a guarantee of feasible-design discovery.

use super::*;

const MAGIC: &[u8] = b"fs-topols/projected-multiload/checkpoint/2\n";

fn validate_reduction(value: f64) -> Result<(), CutFemError> {
    if !(value.is_finite() && (0.0..1.0).contains(&value)) {
        return Err(invalid("stress restoration reduction must be finite and lie in [0,1)"));
    }
    Ok(())
}

// Strict reduction of the worst absolute excess. Reaching the declared bound
// always succeeds; an unchanged/increased stress can never be restoration.
fn reduces_violation(previous: f64, candidate: f64, admitted: f64, reduction: f64) -> bool {
    previous.is_finite() && candidate.is_finite() && admitted.is_finite()
        && previous > admitted && candidate >= 0.0 && admitted > 0.0
        && reduction.is_finite() && (0.0..1.0).contains(&reduction)
        && (candidate <= admitted
            || candidate - admitted < (previous - admitted) * (1.0 - reduction))
}

impl MultiLoadProjectedOptimizer {
    /// Admit an area-feasible but possibly overstressed baseline explicitly.
    ///
    /// Before feasibility, accepted updates strictly reduce the worst sampled
    /// stress excess (by more than `min_relative_reduction` while still above
    /// the bound); compliance may increase and is NOT an acceptance criterion.
    /// Once feasible, the ordinary stress-preserving compliance rule resumes.
    /// The original total update and case-solve limits cover BOTH phases.
    ///
    /// Uses cached baseline displacements, including zero-objective-weight cases.
    /// No stress gradient or new material/allowable is inferred. Configuration
    /// is immutable after installation. The baseline remains explicitly the
    /// original AREA-feasible design, not a fabricated stress-feasible design.
    pub fn with_stress_restoration(
        mut self, limit: SampledStressLimit, min_relative_reduction: f64,
    ) -> Result<Self, CutFemError> {
        validate_reduction(min_relative_reduction)?;
        let limit = SampledStressLimit::new(limit.max_von_mises, limit.absolute_tolerance)?;
        if self.stress_limit.is_some() || self.next_iteration != 0
            || self.solves_started != self.kernel.load_cases.len()
        {
            return Err(invalid("stress restoration must be installed once before candidate solves or updates"));
        }
        let measured = match self.sample_stress_controlled(&self.current, |_, _| {
            ControlFlow::<Infallible>::Continue(())
        })? {
            ControlFlow::Continue(value) => value,
            ControlFlow::Break(never) => match never {},
        };
        self.stress_limit = Some(limit);
        self.restoration_reduction = Some(min_relative_reduction);
        self.baseline_stress = Some(measured.clone());
        self.current_stress = Some(measured);
        Ok(self)
    }

    /// Immutable opt-in restoration policy; absent for the original strict mode.
    #[must_use]
    pub const fn stress_restoration_reduction(&self) -> Option<f64> { self.restoration_reduction }

    /// Accepted stress-restoration updates; not compliance-improvement updates.
    #[must_use]
    pub const fn restoration_updates(&self) -> usize { self.restoration_updates }

    /// Whether the next update must restore stress rather than optimize compliance.
    /// Missing required evidence cannot be mistaken for a feasible state.
    #[must_use]
    pub fn is_restoring_stress(&self) -> bool {
        self.restoration_reduction.is_some()
            && stress::require_feasible(self.stress_limit, self.current_stress.as_ref()).is_err()
    }

    pub(super) fn require_candidate(
        &self, state: &MultiLoadProjectedState, measured: Option<&RobustSampledStressEvaluation>,
    ) -> Result<(), CutFemError> {
        if self.is_restoring_stress() {
            let previous = self.current_stress.as_ref()
                .ok_or_else(|| invalid("restoration requires complete current stress"))?;
            let candidate = measured
                .ok_or_else(|| invalid("restoration requires complete candidate stress"))?;
            let bound = self.stress_limit
                .ok_or_else(|| invalid("restoration requires an immutable stress limit"))?;
            let reduction = self.restoration_reduction
                .ok_or_else(|| invalid("restoration requires an immutable reduction policy"))?;
            if !reduces_violation(previous.worst_sampled_von_mises,
                candidate.worst_sampled_von_mises, bound.admitted_max(), reduction)
            {
                return Err(invalid("insufficient complete-family sampled stress violation reduction"));
            }
        } else {
            stress::require_feasible(self.stress_limit, measured)?;
            let limit = self.current.objective * (1.0 - self.controls.min_relative_improvement);
            if !(state.objective < limit) {
                return Err(invalid("insufficient same-material aggregate decrease"));
            }
        }
        Ok(())
    }

    /// Versioned restart bytes. Original strict/no-stress studies retain their
    /// exact v1 bytes. Restoration studies use v2, carrying the immutable
    /// reduction and accepted restoration count before the retained v1 layout.
    #[must_use]
    pub fn checkpoint_bytes(&self) -> Vec<u8> {
        let payload = self.checkpoint_v1_bytes();
        let Some(reduction) = self.restoration_reduction else { return payload };
        let mut bytes = MAGIC.to_vec();
        bytes.extend_from_slice(&reduction.to_bits().to_le_bytes());
        bytes.extend_from_slice(&(self.restoration_updates as u64).to_le_bytes());
        bytes.extend_from_slice(&payload);
        bytes
    }

    // This checks consistency with both independently re-solved endpoints. It
    // is not authentication of every intermediate step or historical budget use.
    pub(super) fn check_restoration_history(&self) -> Result<(), CutFemError> {
        let Some(_) = self.restoration_reduction else { return Ok(()) };
        let bound = self.stress_limit.ok_or_else(|| invalid("restoration checkpoint lacks its stress limit"))?;
        let baseline = self.baseline_stress.as_ref().ok_or_else(|| invalid("restoration checkpoint lacks baseline stress"))?;
        let current = self.current_stress.as_ref().ok_or_else(|| invalid("restoration checkpoint lacks current stress"))?;
        let initial_infeasible = baseline.worst_sampled_von_mises > bound.admitted_max();
        if (!initial_infeasible && self.restoration_updates != 0)
            || (initial_infeasible && self.next_iteration > 0 && self.restoration_updates == 0)
            || (self.is_restoring_stress() && self.restoration_updates != self.next_iteration)
            || (initial_infeasible && self.next_iteration > 0
                && !(current.worst_sampled_von_mises < baseline.worst_sampled_von_mises))
            || (!initial_infeasible && self.next_iteration > 0
                && !(self.current.objective < self.baseline.objective))
        {
            return Err(invalid("inconsistent restoration checkpoint phase, count or endpoint progress"));
        }
        Ok(())
    }
}

pub(super) fn checkpoint_payload(bytes: &[u8]) -> Result<(&[u8], Option<f64>, usize), CutFemError> {
    if !bytes.starts_with(MAGIC) { return Ok((bytes, None, 0)) }
    let extra = bytes.get(MAGIC.len()..MAGIC.len() + 16)
        .ok_or_else(|| invalid("truncated stress restoration checkpoint"))?;
    let mut word = [0; 8];
    word.copy_from_slice(&extra[..8]);
    let reduction = f64::from_bits(u64::from_le_bytes(word));
    validate_reduction(reduction)?;
    word.copy_from_slice(&extra[8..]);
    let updates = usize::try_from(u64::from_le_bytes(word))
        .map_err(|_| invalid("restoration count exceeds platform range"))?;
    if updates > 10_000 { return Err(invalid("restoration count exceeds the study update cap")); }
    Ok((&bytes[MAGIC.len() + 16..], Some(reduction), updates))
}

#[cfg(test)]
mod tests;

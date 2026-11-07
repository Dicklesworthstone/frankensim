//! Independent local time-error scales for actual finite-rate phase stocks.
//!
//! Water activity is an aggregate dissolved-particle response. It cannot in
//! general resolve each component or a small remaining crystal stock. These
//! optional scales compare the original kinetic owner's dissolved AND solid
//! amounts without replacing its native first-order volume/phase evolution.

use super::*;
use crate::moisture::solution::{FiniteMoistureSolution, SolutionMoistureError};

/// Per-component local temporal accuracy for the actual kinetic phase stocks.
/// Both dissolved and crystalline amounts receive the same absolute allowance.
/// Relative scales use changes from the last accepted state, not the possibly
/// much larger total inventory. This does not control a source/model error.
#[derive(Clone, Copy, Debug)]
pub struct KineticPhaseAccuracy {
    pub absolute_amount_tolerance_mol: f64,
    /// Multiplies max(|coarse-start|, |fine-start|), in [0, 1).
    pub relative_tolerance: f64,
}

/// The existing water/transport controls plus optional kinetic phase accuracy.
/// Existing calls passing [`MoistureTimeControls`] retain their original
/// behavior through `From`; explicit phase control requires an original kinetic
/// reservoir and is admitted before any proposal, including an empty duration.
#[derive(Clone, Copy, Debug)]
pub struct SolutionMoistureTimeControls {
    pub transport: MoistureTimeControls,
    pub kinetic_phases: Option<KineticPhaseAccuracy>,
}

impl From<MoistureTimeControls> for SolutionMoistureTimeControls {
    fn from(transport: MoistureTimeControls) -> Self {
        Self {
            transport,
            kinetic_phases: None,
        }
    }
}

impl MoistureTimeControls {
    /// Resolve actual dissolved/crystal inventories in addition to water,
    /// moisture and activity, using the same bounded full/half controller.
    pub const fn with_kinetic_phase_accuracy(
        self,
        phases: KineticPhaseAccuracy,
    ) -> SolutionMoistureTimeControls {
        SolutionMoistureTimeControls {
            transport: self,
            kinetic_phases: Some(phases),
        }
    }
}

impl SolutionMoistureTimeControls {
    pub(super) fn admit(
        self,
        original: &FiniteMoistureSolution,
        cx: &Cx<'_>,
    ) -> Result<(), MoistureAdvanceError<SolutionMoistureError>> {
        checkpoint(cx)?;
        if let Some(phases) = self.kinetic_phases {
            if !positive(phases.absolute_amount_tolerance_mol)
                || !phases.relative_tolerance.is_finite()
                || !(0.0..1.0).contains(&phases.relative_tolerance)
            {
                return Err(MoistureAdvanceError::Invalid(
                    "positive finite kinetic amount tolerance and relative tolerance in [0,1) required",
                ));
            }
            if original.kinetic_solubility().is_none() {
                return Err(MoistureAdvanceError::Invalid(
                    "kinetic phase accuracy requires the original finite-rate phase owner",
                ));
            }
        }
        Ok(())
    }

    /// Raw full/fine discrepancy; the existing time controller applies its
    /// representable-half correction exactly once after combining all norms.
    pub(super) fn phase_error(
        self,
        old: &FiniteMoistureSolution,
        coarse: &FiniteMoistureSolution,
        fine: &FiniteMoistureSolution,
        cx: &Cx<'_>,
    ) -> Result<f64, MoistureAdvanceError<SolutionMoistureError>> {
        let Some(accuracy) = self.kinetic_phases else {
            return Ok(0.0);
        };
        checkpoint(cx)?;
        let (Some(old), Some(coarse), Some(fine)) = (
            old.kinetic_solubility(),
            coarse.kinetic_solubility(),
            fine.kinetic_solubility(),
        ) else {
            return Err(MoistureAdvanceError::Invalid(
                "native solution proposal changed its kinetic phase owner",
            ));
        };
        let count = old.last_steps().len();
        if count == 0
            || coarse.last_steps().len() != count
            || fine.last_steps().len() != count
            || old.models().len() != count
            || coarse.models().len() != count
            || fine.models().len() != count
        {
            return Err(MoistureAdvanceError::Invalid(
                "native kinetic component layout changed during temporal comparison",
            ));
        }
        let mut ratio = 0.0_f64;
        for i in 0..count {
            checkpoint(cx)?;
            if old.models()[i] != coarse.models()[i] || old.models()[i] != fine.models()[i] {
                return Err(MoistureAdvanceError::Invalid(
                    "native kinetic component source changed during temporal comparison",
                ));
            }
            let a = old.last_steps()[i].state();
            let b = coarse.last_steps()[i].state();
            let c = fine.last_steps()[i].state();
            // Read both native stocks. Reconstructing solid = total-dissolved
            // would erase a real crystal tail below one ulp of the total.
            for [before, full, half] in [
                [
                    a.dissolved().value(),
                    b.dissolved().value(),
                    c.dissolved().value(),
                ],
                [a.solid().value(), b.solid().value(), c.solid().value()],
            ] {
                if [before, full, half]
                    .iter()
                    .any(|n| !n.is_finite() || *n < 0.0)
                {
                    return Err(MoistureAdvanceError::Invalid(
                        "unrepresentable native kinetic phase inventory",
                    ));
                }
                let relative =
                    accuracy.relative_tolerance * (full - before).abs().max((half - before).abs());
                let absolute = accuracy.absolute_amount_tolerance_mol;
                let scale = absolute.max(relative);
                let error = ((half - full).abs() / scale) / (absolute / scale + relative / scale);
                ratio = ratio.max(error);
            }
        }
        Ok(ratio)
    }
}

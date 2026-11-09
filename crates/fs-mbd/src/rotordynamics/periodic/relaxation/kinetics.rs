//! Held dissipative clocks on original Maxwell branches, without history reset.

use super::super::{
    MAX_ROTOR_PERIODIC_SEGMENTS, RotorError, RotorFloquetBudget, RotorPeriodicForcing,
    RotorPeriodicSegment, poll,
};
use super::{
    PreparedRotorMaxwellCycle, RotorMaxwellFloquetResponse,
    power::{PreparedRotorMaxwellPowerCycle, prepare_power_cycle},
    prepare,
};
use crate::rotordynamics::relaxation::{MAX_ROTOR_MAXWELL_BRANCHES, MaxwellDampingBranch};

/// One positive relaxation-time ratio per original [segment][branch].
///
/// Each effective time is `tau_reference * ratio`, ALWAYS from the original
/// supplied branch. Rows are never compounded. Changing a dissipative clock
/// changes neither the branch's stiffness/strain identity nor its stress and
/// stored energy. Every memory continues through each switch and period seam;
/// the actual held time enters both motion and its independently integrated heat.
///
/// This low-level schedule infers no temperature law, calibration support,
/// thermal feedback or smooth-time approximation. A material/thermal owner
/// evaluates its original law and supplies each held row. Changing branch k,
/// strain basis or topology requires a separate storage/history-transfer law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotorMaxwellKinetics<'a> {
    /// Finite positive tau/tau_reference ratios in original segment/branch order.
    pub relaxation_time_scales: &'a [Vec<f64>],
}

impl RotorMaxwellKinetics<'_> {
    /// Validate complete layout and positive representable effective times AND
    /// inverse times without allocation or preparing a flow. The ordinary
    /// physical owner still admits M/G, branch stiffness and strain vectors.
    pub fn validate(
        &self,
        segment_count: usize,
        branches: &[MaxwellDampingBranch],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), RotorError> {
        poll(&mut cancelled)?;
        if segment_count == 0
            || segment_count > MAX_ROTOR_PERIODIC_SEGMENTS
            || branches.len() > MAX_ROTOR_MAXWELL_BRANCHES
            || self.relaxation_time_scales.len() != segment_count
        {
            return Err(RotorError::Invalid(
                "complete bounded Maxwell kinetics segment/branch layout required",
            ));
        }
        for row in self.relaxation_time_scales {
            poll(&mut cancelled)?;
            if row.len() != branches.len() {
                return Err(RotorError::Invalid(
                    "one time ratio per original Maxwell branch required",
                ));
            }
            for (branch, &scale) in branches.iter().zip(row) {
                poll(&mut cancelled)?;
                relaxation_time(branch.relaxation_time_s, Some(scale))?;
            }
        }
        Ok(())
    }

    /// Full homogeneous Floquet map under the authored held material clocks.
    /// Every physical memory multiplier remains; no source oscillator is added.
    pub fn floquet(
        &self,
        spin_rad_s: f64,
        segments: &[RotorPeriodicSegment<'_>],
        branches: &[MaxwellDampingBranch],
        relative_generator_residual_tolerance: f64,
        budget: RotorFloquetBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<RotorMaxwellFloquetResponse, RotorError> {
        Ok(prepare(
            spin_rad_s,
            segments,
            branches,
            Some(*self),
            RotorPeriodicForcing::Homogeneous,
            relative_generator_residual_tolerance,
            budget,
            false,
            cancelled,
        )?
        .report)
    }

    /// Prepare the existing full-state homogeneous, held, harmonic or Fourier
    /// cycle with this schedule. The original force/phase and generator owners
    /// still perform every solve; warm advances/orbits use their unchanged maps.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_cycle(
        &self,
        spin_rad_s: f64,
        segments: &[RotorPeriodicSegment<'_>],
        branches: &[MaxwellDampingBranch],
        forcing: RotorPeriodicForcing<'_>,
        relative_generator_residual_tolerance: f64,
        budget: RotorFloquetBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedRotorMaxwellCycle, RotorError> {
        prepare(
            spin_rad_s,
            segments,
            branches,
            Some(*self),
            forcing,
            relative_generator_residual_tolerance,
            budget,
            true,
            cancelled,
        )
    }

    /// Prepare the same native motion plus original source work, storage and
    /// every branch's heat using its actual interval time. A pure clock switch
    /// has no storage jump or modulation work; ordinary relaxed-K switches keep
    /// their independently observed work. A homogeneous source is valid too.
    /// The schedule is consumed only during preparation, never borrowed by the
    /// resulting cache. Original storage laws remain borrowed and unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_power_cycle<'a>(
        &self,
        spin_rad_s: f64,
        segments: &[RotorPeriodicSegment<'a>],
        branches: &'a [MaxwellDampingBranch],
        forcing: RotorPeriodicForcing<'_>,
        relative_generator_residual_tolerance: f64,
        budget: RotorFloquetBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedRotorMaxwellPowerCycle<'a>, RotorError> {
        prepare_power_cycle(
            spin_rad_s,
            segments,
            branches,
            Some(*self),
            forcing,
            relative_generator_residual_tolerance,
            budget,
            cancelled,
        )
    }
}

// All interval times derive from immutable original branch data. Keeping this
// one owner for admission, generator, independent residual and heat prevents a
// history or energetic discrepancy when a material clock changes.
pub(super) fn relaxation_time(reference_s: f64, scale: Option<f64>) -> Result<f64, RotorError> {
    let Some(scale) = scale else {
        return Ok(reference_s);
    };
    let time = reference_s * scale;
    let rate = time.recip();
    if !reference_s.is_finite()
        || reference_s <= 0.0
        || !scale.is_finite()
        || scale <= 0.0
        || !time.is_finite()
        || time <= 0.0
        || !rate.is_finite()
        || rate <= 0.0
    {
        Err(RotorError::Invalid(
            "Maxwell clock ratios require positive representable times and rates",
        ))
    } else {
        Ok(time)
    }
}

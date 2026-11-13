//! Native curved-shell endpoints share the bounded whole-duration controller.

use super::*;
use crate::moisture::radial::MoistureRadial;

mod finite;
mod solution;

fn snapshot(state: &MoistureRadial) -> Result<Snapshot<'_>, MoistureError> {
    Ok(Snapshot {
        fields: Specimens::One((0, state.field())),
        field_activities: None,
        reservoir_water_kg: None,
        activity: None,
    })
}

/// An L1 water discrepancy retains opposing shell errors: a closed specimen's
/// zero net water change must not hide an inaccurate internal redistribution.
/// Form local differences before multiplying by dry mass, avoiding subtraction
/// of large stored totals. Relative tolerance scales the actual transported
/// water, not the potentially much larger water already bound in the specimen.
pub(super) fn water_distance<E>(
    original: &MoistureRadial,
    coarse: &MoistureRadial,
    fine: &MoistureRadial,
    controls: MoistureTimeControls,
    cx: &Cx<'_>,
) -> Result<f64, MoistureAdvanceError<E>> {
    let mut difference = WaterSum::default();
    let mut coarse_change = WaterSum::default();
    let mut fine_change = WaterSum::default();
    for (i, ((&before, &a), &b)) in original
        .moisture()
        .iter()
        .zip(coarse.moisture())
        .zip(fine.moisture())
        .enumerate()
    {
        if i % 128 == 0 {
            checkpoint(cx)?;
        }
        let mass = original.dry_masses_kg()[i];
        difference.add(mass * (a - b).abs());
        coarse_change.add(mass * (a - before).abs());
        fine_change.add(mass * (b - before).abs());
    }
    let changes = [coarse_change.total(), fine_change.total()];
    let allowance = controls.reservoir_water_absolute_tolerance_kg
        + controls.relative_tolerance * changes[0].max(changes[1]);
    let error = difference.total();
    if changes.iter().any(|change| !change.is_finite())
        || !positive(allowance)
        || !error.is_finite()
        || error < 0.0
    {
        return Err(MoistureAdvanceError::Invalid(
            "unrepresentable radial temporal water comparison",
        ));
    }
    Ok(error / allowance)
}

impl MoistureRadial {
    /// Evolve a complete duration with nonlinear backward Euler on the original
    /// cylindrical or spherical shells. The original endpoint source law,
    /// curved resistance, center condition and surface exchange stay in use.
    ///
    /// Acceptance resolves every shell's moisture and the sum of absolute
    /// shell water discrepancies. For this radial owner,
    /// `reservoir_water_absolute_tolerance_kg` also sets that whole-specimen
    /// water error scale; opposing interior transfers cannot cancel it. Actual
    /// per-shell dry masses weight both the comparison and final water account.
    /// Relative tolerance scales changes from the last accepted state.
    ///
    /// Faces and source temperature stay fixed. The full-versus-two-half error
    /// estimate is local and asymptotic; it is not a global accuracy guarantee.
    /// Only accepted half-step fluxes enter the fixed whole-duration water
    /// budget. Any failure discards the entire proposed duration, including
    /// successful internal prefixes.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive_to(
        &self,
        target_s: f64,
        faces: &[MoistureFace; 2],
        step: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (Self, MoistureAdvanceReport<MoistureStepReport, 2>),
        MoistureAdvanceError<MoistureError>,
    > {
        let (state, report) = evolve_fields_with_method(
            self,
            target_s,
            time,
            TimeMethod::FirstOrder,
            |state, dt| {
                state
                    .advance(dt, faces, step, nonlinear, cx)
                    .map(|(state, report)| (state, OneReport(report)))
            },
            snapshot,
            |old, coarse, fine| water_distance(old, coarse, fine, time, cx),
            cx,
        )?;
        Ok((state, report.summary.map_step(|step| step.0)))
    }

    /// Evolve a complete duration with the original second-order SDIRK2 stages.
    /// The same componentwise moisture and dry-mass-weighted absolute water
    /// comparisons as [`Self::advance_adaptive_to`] apply. The controller uses
    /// cubic local-error scaling and the correct unequal-half-step factor.
    ///
    /// Each attempted endpoint charges two implicit stages against the work
    /// limit, including an early failure; a full comparison costs six stages.
    /// Physical `StageRange` refusals shorten the step without clipping or a
    /// first-order fallback. Only accepted RK-weighted boundary transfers enter
    /// the original-to-final conservation account. A later failure or
    /// cancellation publishes no accepted prefix.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_sdirk2_adaptive_to(
        &self,
        target_s: f64,
        faces: &[MoistureFace; 2],
        step: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (Self, MoistureAdvanceReport<MoistureSdirk2Report, 2>),
        MoistureAdvanceError<MoistureError>,
    > {
        let (state, report) = evolve_fields_with_method(
            self,
            target_s,
            time,
            TimeMethod::Sdirk2,
            |state, dt| {
                state
                    .advance_sdirk2(dt, faces, step, nonlinear, cx)
                    .map(|(state, report)| (state, OneReport(report)))
            },
            snapshot,
            |old, coarse, fine| water_distance(old, coarse, fine, time, cx),
            cx,
        )?;
        Ok((state, report.summary.map_step(|step| step.0)))
    }
}

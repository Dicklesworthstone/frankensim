//! Whole-duration curved diffusion and the original finite aqueous inventory.

use super::*;
use crate::moisture::{
    atmosphere::DryBasisSorption,
    solution::{
        FiniteMoistureSolution, SolutionMoistureError, SolutionMoistureSolve, SolutionMoistureStep,
    },
};

fn solution_snapshot(
    state: &(FiniteMoistureSolution, MoistureRadial),
) -> Result<Snapshot<'_>, SolutionMoistureError> {
    Ok(Snapshot {
        fields: Specimens::One((0, state.1.field())),
        field_activities: None,
        reservoir_water_kg: Some(state.0.solvent_mass_kg()),
        activity: Some(state.0.solution().water_activity()?),
    })
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
impl FiniteMoistureSolution {
    /// Complete radial uptake/desorption and finite-solution feedback with
    /// backward-Euler step doubling. Every proposal solves the original solvent
    /// volume/activity root with the source-owned curved diffusion and faces.
    /// No fixed-activity bath or equal-mass shell approximation is introduced.
    ///
    /// Acceptance checks each content, the actual dry-mass-weighted L1 shell
    /// water discrepancy, solvent mass and solution activity. The water scale
    /// applies to both the shell sum and solvent comparison; pure solvent is
    /// still resolved when its activity remains one. Relative scales refer to
    /// changes since the last accepted state. Only accepted fine exchanges enter
    /// the original-to-final water and transport accounts.
    ///
    /// The original bounded retries and source/phase admission remain in use.
    /// Kinetic phases retain their native first-order endpoint law. Passing
    /// `time.with_kinetic_phase_accuracy(...)` additionally controls each actual
    /// dissolved/crystalline amount, including small remaining crystal stocks.
    /// Temperature and faces are held fixed, and the original equal-pressure
    /// vapor bridge does not imply liquid immersion. Any later failure discards
    /// both endpoints, including successful private prefixes. Local estimates
    /// do not assert a global accuracy bound.
    pub fn advance_radial_adaptive_to(
        &self,
        moisture: &MoistureRadial,
        target_s: f64,
        faces: &[Option<DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: SolutionMoistureSolve,
        time: impl Into<SolutionMoistureTimeControls>,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            MoistureRadial,
            MoistureAdvanceReport<SolutionMoistureStep, 2>,
        ),
        MoistureAdvanceError<SolutionMoistureError>,
    > {
        checkpoint(cx)?;
        let accuracy = time.into();
        accuracy.admit(self, cx)?;
        let time = accuracy.transport;
        let initial = (self.clone(), moisture.clone());
        let ((reservoir, field), report) = evolve_fields_with_method(
            &initial,
            target_s,
            time,
            TimeMethod::FirstOrder,
            |state, dt| {
                let (bath, field, step) = state
                    .0
                    .advance_radial(&state.1, dt, faces, transport, nonlinear, physical, cx)?
                    .into_parts();
                Ok(((bath, field), OneReport(step)))
            },
            solution_snapshot,
            |old, coarse, fine| {
                Ok(water_distance(&old.1, &coarse.1, &fine.1, time, cx)?
                    .max(accuracy.phase_error(&old.0, &coarse.0, &fine.0, cx)?))
            },
            cx,
        )?;
        Ok((reservoir, field, report.summary.map_step(|step| step.0)))
    }

    /// Complete a finite-solution radial duration with coupled SDIRK2 stages.
    /// Each stage solves its own actual bath activity, sorption, surface flux
    /// and source-dependent shell diffusion. The endpoint includes the complete
    /// RK-weighted water exchange; diagnostic internal stages are not debited
    /// again. The original controller uses second-order step doubling and cubic
    /// local-error scaling with the comparisons in [`Self::advance_radial_adaptive_to`].
    ///
    /// Both implicit stages are reserved against the work cap even on failure.
    /// Range/convergence errors shorten proposals within original bounds,
    /// without clipping, expanding source domains or falling back to BE.
    /// Algebraic phase partitioning may be used on its admitted smooth branch;
    /// kinetic phase inventories remain refused by the native SDIRK owner.
    /// No separate phase-component error estimate or phase heat is supplied.
    /// Final solvent and shell inventories plus accepted face exchanges retain
    /// one fixed whole-duration water tolerance before atomic publication.
    pub fn advance_radial_sdirk2_adaptive_to(
        &self,
        moisture: &MoistureRadial,
        target_s: f64,
        faces: &[Option<DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: SolutionMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            MoistureRadial,
            MoistureAdvanceReport<SolutionMoistureStep<MoistureSdirk2Report>, 2>,
        ),
        MoistureAdvanceError<SolutionMoistureError>,
    > {
        checkpoint(cx)?;
        let initial = (self.clone(), moisture.clone());
        let ((reservoir, field), report) = evolve_fields_with_method(
            &initial,
            target_s,
            time,
            TimeMethod::Sdirk2,
            |state, dt| {
                let (bath, field, step) = state
                    .0
                    .advance_radial_sdirk2(&state.1, dt, faces, transport, nonlinear, physical, cx)?
                    .into_parts();
                Ok(((bath, field), OneReport(step)))
            },
            solution_snapshot,
            |old, coarse, fine| water_distance(&old.1, &coarse.1, &fine.1, time, cx),
            cx,
        )?;
        Ok((reservoir, field, report.summary.map_step(|step| step.0)))
    }
}

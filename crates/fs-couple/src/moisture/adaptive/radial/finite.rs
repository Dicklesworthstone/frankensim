//! Whole-duration radial diffusion and one original finite vapor inventory.

use super::*;
use crate::moisture::{
    adaptive::finite::check_vapor,
    atmosphere::{
        AtmosphereMoistureError, AtmosphereMoistureSolve, AtmosphereMoistureStep, DryBasisSorption,
        FiniteMoistureAtmosphere,
    },
};

fn vapor_snapshot(
    state: &(FiniteMoistureAtmosphere, MoistureRadial),
) -> Result<Snapshot<'_>, AtmosphereMoistureError> {
    Ok(Snapshot {
        fields: Specimens::One((0, state.1.field())),
        field_activities: None,
        reservoir_water_kg: Some(state.0.vapor_mass_kg()?),
        activity: Some(state.0.relative_humidity()?),
    })
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
impl FiniteMoistureAtmosphere {
    /// Evolve one finite atmosphere and its original cylinder or sphere to a
    /// common final time using nonlinear backward Euler. Every probe solves the
    /// original coupled endpoint humidity, sorption and curved diffusion laws.
    /// No fixed-humidity bath substitutes for the evolving gas inventory.
    ///
    /// Temporal acceptance checks each shell's content, the dry-mass-weighted
    /// L1 water discrepancy, actual vapor mass and RH. The water tolerance sets
    /// both the shell sum and vapor comparison scales. Only accepted fine
    /// endpoints contribute boundary exchange; rejected trials debit nothing.
    ///
    /// Whole-duration checks retain separate water, carrier and total-mass
    /// accounts with one fixed budget. Temperature and gas volume stay fixed;
    /// gas energy change is not the total sorption heat. A later failure or
    /// cancellation returns neither proposed inventory, including any accepted
    /// private prefix. Local error estimates are not global accuracy bounds.
    pub fn advance_radial_adaptive_to(
        &self,
        moisture: &MoistureRadial,
        target_s: f64,
        faces: &[Option<DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: AtmosphereMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            MoistureRadial,
            MoistureAdvanceReport<AtmosphereMoistureStep, 2>,
        ),
        MoistureAdvanceError<AtmosphereMoistureError>,
    > {
        checkpoint(cx)?;
        let initial = (self.clone(), moisture.clone());
        let ((atmosphere, field), report) = evolve_fields_with_method(
            &initial,
            target_s,
            time,
            TimeMethod::FirstOrder,
            |state, dt| {
                let (air, field, step) = state
                    .0
                    .advance_radial(&state.1, dt, faces, transport, nonlinear, physical, cx)?
                    .into_parts();
                Ok(((air, field), OneReport(step)))
            },
            vapor_snapshot,
            |old, coarse, fine| water_distance(&old.1, &coarse.1, &fine.1, time, cx),
            cx,
        )?;
        let mut report = report.summary.map_step(|step| step.0);
        check_vapor(
            self,
            &atmosphere,
            report.water_change_kg,
            time.water_balance_tolerance_kg,
            &mut report,
            cx,
        )?;
        Ok((atmosphere, field, report))
    }

    /// Evolve one finite atmosphere and its original radial field with the
    /// coupled SDIRK2 endpoint at every full and half proposal. Each native stage
    /// solves its own actual gas humidity and source-dependent shell diffusion;
    /// the controller never substitutes prescribed-humidity transport.
    ///
    /// The same moisture, local water, vapor and RH comparisons as
    /// [`Self::advance_radial_adaptive_to`] use second-order step doubling and
    /// cubic local-error scaling. Each attempted endpoint reserves two implicit
    /// stages, even on early failure. Only accepted RK-weighted face exchanges
    /// are accumulated, so diagnostic internal stages are not debited twice.
    ///
    /// Physical range and convergence refusals shorten proposals within the
    /// original limits, without clipping, widening sources or first-order
    /// fallback. Both physical endpoints remain private until final water,
    /// carrier and total-mass checks and cancellation have succeeded.
    pub fn advance_radial_sdirk2_adaptive_to(
        &self,
        moisture: &MoistureRadial,
        target_s: f64,
        faces: &[Option<DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: AtmosphereMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            MoistureRadial,
            MoistureAdvanceReport<AtmosphereMoistureStep<MoistureSdirk2Report>, 2>,
        ),
        MoistureAdvanceError<AtmosphereMoistureError>,
    > {
        checkpoint(cx)?;
        let initial = (self.clone(), moisture.clone());
        let ((atmosphere, field), report) = evolve_fields_with_method(
            &initial,
            target_s,
            time,
            TimeMethod::Sdirk2,
            |state, dt| {
                let (air, field, step) = state
                    .0
                    .advance_radial_sdirk2(&state.1, dt, faces, transport, nonlinear, physical, cx)?
                    .into_parts();
                Ok(((air, field), OneReport(step)))
            },
            vapor_snapshot,
            |old, coarse, fine| water_distance(&old.1, &coarse.1, &fine.1, time, cx),
            cx,
        )?;
        let mut report = report.summary.map_step(|step| step.0);
        check_vapor(
            self,
            &atmosphere,
            report.water_change_kg,
            time.water_balance_tolerance_kg,
            &mut report,
            cx,
        )?;
        Ok((atmosphere, field, report))
    }
}

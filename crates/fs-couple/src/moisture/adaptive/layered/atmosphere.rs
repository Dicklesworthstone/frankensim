//! Native two-stage finite-air diffusion under the original duration controller.

use super::radial::{observe_radial, radial_water_distance};
use super::*;
use crate::moisture::radial::layered::MoistureLayeredRadial;

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
impl FiniteMoistureAtmosphere {
    /// Complete heterogeneous slab diffusion and one actual finite atmosphere
    /// with native coupled SDIRK2 endpoints. Each stage solves gas humidity and
    /// every original nonlinear material/interface law together.
    ///
    /// Local moisture/activity, actual vapor mass/RH and the L1 dry-mass-weighted
    /// water discrepancy enter the original second-order error comparison.
    /// Both native stages spend the work cap, including rejected probes. Only
    /// accepted fine endpoint transfers enter the material and exterior accounts;
    /// their first-stage diagnostics are never debited a second time.
    ///
    /// Original-to-final water, total gas/solid mass and each carrier are checked
    /// over the complete duration against one fixed tolerance. A late refusal or
    /// cancellation publishes no accepted prefix. Temperature, gas free volume,
    /// reference geometry and exposure stay fixed. Local temporal estimates do
    /// not certify global accuracy or add sorption heat or condensation.
    pub fn advance_layered_sdirk2_adaptive_to(
        &self,
        moisture: &MoistureLayeredSlab,
        target_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: AtmosphereMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            MoistureLayeredSlab,
            LayeredMoistureAdvanceReport<AtmosphereMoistureStep<MoistureLayeredSdirk2Report>>,
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
                    .advance_layered_sdirk2(
                        &state.1, dt, exposed, transport, nonlinear, physical, cx,
                    )?
                    .into_parts();
                Ok(((air, field), LayeredStep(step)))
            },
            |state| {
                let mut snapshot = observe(&state.1);
                snapshot.reservoir_water_kg = Some(state.0.vapor_mass_kg()?);
                snapshot.activity = Some(state.0.relative_humidity()?);
                Ok(snapshot)
            },
            |old, coarse, fine| {
                water_distance(
                    old.1
                        .layers()
                        .iter()
                        .zip(coarse.1.layers())
                        .zip(fine.1.layers())
                        .map(|((a, b), c)| (a.slab().field(), b.slab().field(), c.slab().field())),
                    time,
                    cx,
                )
            },
            cx,
        )?;
        let mut report = finish(report, time.water_balance_tolerance_kg, cx)?;
        finite::check_vapor(
            self,
            &atmosphere,
            report.summary.water_change_kg,
            time.water_balance_tolerance_kg,
            &mut report.summary,
            cx,
        )?;
        Ok((atmosphere, field, report))
    }

    /// Complete coated cylinder/sphere diffusion with the same finite-air
    /// SDIRK stages and whole-duration checks as
    /// [`Self::advance_layered_sdirk2_adaptive_to`]. Every shell keeps its actual
    /// unequal dry mass and curved resistance. Material joins have paired
    /// transfers; only the exposed inner/outer physical surfaces debit gas.
    ///
    /// Local water errors cannot cancel across layers, even with sealed
    /// exteriors. A zero-area center cannot be exposed. Any final gas or material
    /// account refusal discards both complete candidates.
    pub fn advance_layered_radial_sdirk2_adaptive_to(
        &self,
        moisture: &MoistureLayeredRadial,
        target_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: AtmosphereMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            MoistureLayeredRadial,
            LayeredMoistureAdvanceReport<AtmosphereMoistureStep<MoistureLayeredSdirk2Report>>,
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
                    .advance_layered_radial_sdirk2(
                        &state.1, dt, exposed, transport, nonlinear, physical, cx,
                    )?
                    .into_parts();
                Ok(((air, field), LayeredStep(step)))
            },
            |state| {
                let mut snapshot = observe_radial(&state.1);
                snapshot.reservoir_water_kg = Some(state.0.vapor_mass_kg()?);
                snapshot.activity = Some(state.0.relative_humidity()?);
                Ok(snapshot)
            },
            |old, coarse, fine| radial_water_distance(&old.1, &coarse.1, &fine.1, time, cx),
            cx,
        )?;
        let mut report = finish(report, time.water_balance_tolerance_kg, cx)?;
        finite::check_vapor(
            self,
            &atmosphere,
            report.summary.water_change_kg,
            time.water_balance_tolerance_kg,
            &mut report.summary,
            cx,
        )?;
        Ok((atmosphere, field, report))
    }
}

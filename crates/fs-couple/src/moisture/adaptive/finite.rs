//! Typed finite-reservoir bindings to the one moisture temporal controller.

use super::*;
use crate::moisture::{
    atmosphere::{
        AtmosphereMoistureError, AtmosphereMoistureSolve, AtmosphereMoistureStep,
        CartesianAtmosphereMoistureStep, CartesianAtmosphereSdirk2Step, DryBasisSorption,
        FiniteMoistureAtmosphere, SlabAtmosphereSdirk2Step,
    },
    solution::{
        CartesianSolutionMoistureStep, CartesianSolutionSdirk2Step, FiniteMoistureSolution,
        SlabSolutionSdirk2Step, SolutionMoistureError, SolutionMoistureSolve, SolutionMoistureStep,
    },
};

impl Retryable for AtmosphereMoistureError {
    fn retryable(&self) -> bool {
        match self {
            Self::Moisture(error) => error.retryable(),
            Self::NotConverged { .. } => true,
            _ => false,
        }
    }
}
impl Retryable for SolutionMoistureError {
    fn retryable(&self) -> bool {
        match self {
            Self::Moisture(error) => error.retryable(),
            Self::NotConverged { .. } => true,
            _ => false,
        }
    }
}
impl<R: WaterStep<N>, const N: usize> WaterStep<N> for AtmosphereMoistureStep<R, N> {
    fn boundary_water(&self) -> [f64; N] {
        self.transport.boundary_water()
    }
}
impl<R: WaterStep<N>, const N: usize> WaterStep<N> for SolutionMoistureStep<R> {
    fn boundary_water(&self) -> [f64; N] {
        self.transport.boundary_water()
    }
}

fn vapor_snapshot<F: ReservoirField<N>, const N: usize>(
    state: &(FiniteMoistureAtmosphere, F),
) -> Result<Snapshot<'_>, AtmosphereMoistureError> {
    Ok(Snapshot {
        fields: Specimens::One((0, state.1.inventory())),
        field_activities: None,
        reservoir_water_kg: Some(state.0.vapor_mass_kg()?),
        activity: Some(state.0.relative_humidity()?),
    })
}
fn solution_snapshot<F: ReservoirField<N>, const N: usize>(
    state: &(FiniteMoistureSolution, F),
) -> Result<Snapshot<'_>, SolutionMoistureError> {
    Ok(Snapshot {
        fields: Specimens::One((0, state.1.inventory())),
        field_activities: None,
        reservoir_water_kg: Some(state.0.solvent_mass_kg()),
        activity: Some(state.0.solution().water_activity()?),
    })
}
pub(super) fn check_vapor<R, const N: usize>(
    original: &FiniteMoistureAtmosphere,
    next: &FiniteMoistureAtmosphere,
    water_change: f64,
    tolerance: f64,
    report: &mut MoistureAdvanceReport<R, N>,
    cx: &Cx<'_>,
) -> Result<(), MoistureAdvanceError<AtmosphereMoistureError>> {
    let residuals = original
        .exchange_residuals(next, water_change)
        .map_err(MoistureAdvanceError::FinalInventory)?;
    for residual in residuals {
        if !residual.is_finite() || residual.abs() > tolerance {
            return Err(MoistureAdvanceError::FinalInventory(
                AtmosphereMoistureError::MassBalance {
                    residual_kg: residual,
                    tolerance_kg: tolerance,
                },
            ));
        }
    }
    report.maximum_carrier_mass_residual_kg = Some(residuals[1]);
    report.total_mass_residual_kg = Some(residuals[2]);
    checkpoint(cx)
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
impl FiniteMoistureAtmosphere {
    /// Complete a positive physical duration with one evolving finite vapor
    /// inventory and the original slab. Faces/temperature stay fixed; source
    /// RH and optional surface emission refresh at every native endpoint.
    /// Temporal error checks content, RH and vapor mass independently. Final
    /// water, total gas/solid mass and each carrier retain a fixed whole-window
    /// budget; no failed/rejected probe is debited or returned.
    pub fn advance_adaptive_to(
        &self,
        moisture: &MoistureSlab,
        target_s: f64,
        faces: &[Option<DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        physical: AtmosphereMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            MoistureSlab,
            MoistureAdvanceReport<AtmosphereMoistureStep, 2>,
        ),
        MoistureAdvanceError<AtmosphereMoistureError>,
    > {
        checkpoint(cx)?;
        let initial = (self.clone(), moisture.clone());
        let ((atmosphere, field), mut report) = evolve(
            &initial,
            target_s,
            time,
            |state, dt| {
                let (air, field, step) = state
                    .0
                    .advance(&state.1, dt, faces, transport, physical, cx)?
                    .into_parts();
                Ok(((air, field), step))
            },
            vapor_snapshot::<MoistureSlab, 2>,
            cx,
        )?;
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

    /// Complete a duration with the native slab's second-order diffusion and
    /// actual finite vapor inventory at every implicit stage.
    ///
    /// The existing temporal controller compares one full interval with two
    /// physical half intervals, uses the second-order error factor (one third
    /// for equal halves) and the cubic growth rule. Moisture in every cell, RH
    /// and vapor mass have independent error scales; only accepted half-step transfers
    /// enter the fixed whole-duration water/carrier/total-mass checks.
    ///
    /// Both implicit stages are charged for every attempted native call.
    /// A retry never publishes its coarse state or diagnostic first-stage gas.
    /// Source refusal, cancellation or a late work limit discards the complete
    /// proposed duration. This local asymptotic estimate is not a global error
    /// bound; source support and the isothermal physical scope are unchanged.
    pub fn advance_sdirk2_adaptive_to(
        &self,
        moisture: &MoistureSlab,
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
            MoistureSlab,
            MoistureAdvanceReport<SlabAtmosphereSdirk2Step, 2>,
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
                    .advance_sdirk2(&state.1, dt, faces, transport, nonlinear, physical, cx)?
                    .into_parts();
                Ok(((air, field), OneReport(step)))
            },
            vapor_snapshot::<MoistureSlab, 2>,
            |_, _, _| Ok(0.0),
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

    /// The same finite-vapor temporal control for the complete Cartesian field,
    /// with six original faces and endpoint-dependent principal diffusivity.
    /// Both physical inventories are returned only after the entire duration.
    pub fn advance_cartesian_adaptive_to(
        &self,
        moisture: &MoistureCartesian,
        target_s: f64,
        faces: &[Option<DryBasisSorption>; 6],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: AtmosphereMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            MoistureCartesian,
            MoistureAdvanceReport<CartesianAtmosphereMoistureStep, 6>,
        ),
        MoistureAdvanceError<AtmosphereMoistureError>,
    > {
        checkpoint(cx)?;
        let initial = (self.clone(), moisture.clone());
        let ((atmosphere, field), mut report) = evolve(
            &initial,
            target_s,
            time,
            |state, dt| {
                let (air, field, step) = state
                    .0
                    .advance_cartesian(&state.1, dt, faces, transport, nonlinear, physical, cx)?
                    .into_parts();
                Ok(((air, field), step))
            },
            vapor_snapshot::<MoistureCartesian, 6>,
            cx,
        )?;
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

    /// Complete a duration with second-order diffusion and one evolving finite
    /// vapor inventory. Every native SDIRK stage solves its own source-backed
    /// RH, sorption and surface transfer together with the Cartesian field.
    ///
    /// The local error estimate checks content, vapor mass and RH separately,
    /// with the second-order step-doubling factor and cubic growth rule. Both
    /// implicit stages are charged even on rejected probes. Only accepted
    /// physical half intervals contribute their RK-weighted face transfers;
    /// diagnostic internal stages are never debited a second time.
    ///
    /// Native source/bracket admissibility and the isothermal model's scope
    /// remain unchanged. A failed duration returns neither physical endpoint,
    /// including when a work limit follows successful private intervals.
    pub fn advance_cartesian_sdirk2_adaptive_to(
        &self,
        moisture: &MoistureCartesian,
        target_s: f64,
        faces: &[Option<DryBasisSorption>; 6],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: AtmosphereMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            MoistureCartesian,
            MoistureAdvanceReport<CartesianAtmosphereSdirk2Step, 6>,
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
                    .advance_cartesian_sdirk2(
                        &state.1, dt, faces, transport, nonlinear, physical, cx,
                    )?
                    .into_parts();
                Ok(((air, field), OneReport(step)))
            },
            vapor_snapshot::<MoistureCartesian, 6>,
            |_, _, _| Ok(0.0),
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

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
impl FiniteMoistureSolution {
    /// Complete a duration through the original equal-pressure solution/vapor
    /// activity bridge. Nonvolatile solute amounts stay fixed; the same finite
    /// solvent inventory evolves along both fine half steps. Local errors check
    /// solvent mass even for pure solvent whose activity remains exactly one.
    /// Pass `time.with_kinetic_phase_accuracy(...)` to additionally resolve each
    /// actual dissolved/crystal stock on the original finite-rate phase owner.
    pub fn advance_adaptive_to(
        &self,
        moisture: &MoistureSlab,
        target_s: f64,
        faces: &[Option<DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        physical: SolutionMoistureSolve,
        time: impl Into<SolutionMoistureTimeControls>,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            MoistureSlab,
            MoistureAdvanceReport<SolutionMoistureStep, 2>,
        ),
        MoistureAdvanceError<SolutionMoistureError>,
    > {
        checkpoint(cx)?;
        let accuracy = time.into();
        accuracy.admit(self, cx)?;
        let time = accuracy.transport;
        let initial = (self.clone(), moisture.clone());
        let ((reservoir, field), report) = evolve_with_error(
            &initial,
            target_s,
            time,
            |state, dt| {
                let (bath, field, step) = state
                    .0
                    .advance(&state.1, dt, faces, transport, physical, cx)?
                    .into_parts();
                Ok(((bath, field), step))
            },
            solution_snapshot::<MoistureSlab, 2>,
            |old, coarse, fine| accuracy.phase_error(&old.0, &coarse.0, &fine.0, cx),
            cx,
        )?;
        Ok((reservoir, field, report))
    }

    /// Complete a physical duration with native second-order slab diffusion
    /// and the same evolving finite solution inventory at every implicit stage.
    ///
    /// The existing controller checks each cell, solvent mass and activity
    /// separately. Pure solvent still has a finite water inventory even though
    /// its activity is constant. The second-order step-doubling estimator and
    /// cubic growth rule admit only the two physical half-step endpoints; their
    /// weighted transfers alone enter the fixed whole-duration water account.
    ///
    /// Each proposal charges both implicit stages. Source refusal, cancellation
    /// or a late work limit discards all private progress. Original algebraic
    /// phase/source admission is unchanged; finite-rate phase stocks remain
    /// explicitly refused. Local smooth-branch accuracy is not a global error
    /// certificate or a guarantee of order across phase transitions.
    pub fn advance_sdirk2_adaptive_to(
        &self,
        moisture: &MoistureSlab,
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
            MoistureSlab,
            MoistureAdvanceReport<SlabSolutionSdirk2Step, 2>,
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
                    .advance_sdirk2(&state.1, dt, faces, transport, nonlinear, physical, cx)?
                    .into_parts();
                Ok(((bath, field), OneReport(step)))
            },
            solution_snapshot::<MoistureSlab, 2>,
            |_, _, _| Ok(0.0),
            cx,
        )?;
        Ok((reservoir, field, report.summary.map_step(|step| step.0)))
    }

    /// Complete Cartesian uptake/desorption and finite-solution concentration
    /// with one temporal controller and the original six-face diffusion owner.
    /// The original activity/volume bracket must support every accepted state.
    pub fn advance_cartesian_adaptive_to(
        &self,
        moisture: &MoistureCartesian,
        target_s: f64,
        faces: &[Option<DryBasisSorption>; 6],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: SolutionMoistureSolve,
        time: impl Into<SolutionMoistureTimeControls>,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            MoistureCartesian,
            MoistureAdvanceReport<CartesianSolutionMoistureStep, 6>,
        ),
        MoistureAdvanceError<SolutionMoistureError>,
    > {
        checkpoint(cx)?;
        let accuracy = time.into();
        accuracy.admit(self, cx)?;
        let time = accuracy.transport;
        let initial = (self.clone(), moisture.clone());
        let ((reservoir, field), report) = evolve_with_error(
            &initial,
            target_s,
            time,
            |state, dt| {
                let (bath, field, step) = state
                    .0
                    .advance_cartesian(&state.1, dt, faces, transport, nonlinear, physical, cx)?
                    .into_parts();
                Ok(((bath, field), step))
            },
            solution_snapshot::<MoistureCartesian, 6>,
            |old, coarse, fine| accuracy.phase_error(&old.0, &coarse.0, &fine.0, cx),
            cx,
        )?;
        Ok((reservoir, field, report))
    }

    /// Complete a duration with second-order diffusion and finite solution
    /// feedback at every implicit stage. Content, solvent mass and water
    /// activity have independent local error scales; the accepted field and
    /// finite bath retain one fixed original-to-final water budget.
    ///
    /// The existing controller uses second-order step-doubling and charges
    /// both implicit stages, including rejected probes. Only accepted physical
    /// half intervals enter the cumulative face ledger. Native source support,
    /// phase admissibility and the equal-pressure, prescribed-temperature
    /// model scope remain unchanged. Kinetic phase inventories are refused by
    /// the native two-stage owner. No partial duration is published on failure.
    pub fn advance_cartesian_sdirk2_adaptive_to(
        &self,
        moisture: &MoistureCartesian,
        target_s: f64,
        faces: &[Option<DryBasisSorption>; 6],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: SolutionMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            MoistureCartesian,
            MoistureAdvanceReport<CartesianSolutionSdirk2Step, 6>,
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
                    .advance_cartesian_sdirk2(
                        &state.1, dt, faces, transport, nonlinear, physical, cx,
                    )?
                    .into_parts();
                Ok(((bath, field), OneReport(step)))
            },
            solution_snapshot::<MoistureCartesian, 6>,
            |_, _, _| Ok(0.0),
            cx,
        )?;
        Ok((reservoir, field, report.summary.map_step(|step| step.0)))
    }
}

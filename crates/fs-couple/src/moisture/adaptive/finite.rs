//! Typed finite-reservoir bindings to the one moisture temporal controller.

use super::*;
use crate::moisture::{
    atmosphere::{
        AtmosphereMoistureError, AtmosphereMoistureSolve, AtmosphereMoistureStep,
        CartesianAtmosphereMoistureStep, DryBasisSorption, FiniteMoistureAtmosphere,
    },
    solution::{
        CartesianSolutionMoistureStep, FiniteMoistureSolution, SolutionMoistureError,
        SolutionMoistureSolve, SolutionMoistureStep,
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
        field: state.1.inventory(),
        reservoir_water_kg: Some(state.0.vapor_mass_kg()?),
        activity: Some(state.0.relative_humidity()?),
    })
}
fn solution_snapshot<F: ReservoirField<N>, const N: usize>(
    state: &(FiniteMoistureSolution, F),
) -> Result<Snapshot<'_>, SolutionMoistureError> {
    Ok(Snapshot {
        field: state.1.inventory(),
        reservoir_water_kg: Some(state.0.solvent_mass_kg()),
        activity: Some(state.0.solution().water_activity()?),
    })
}
fn check_vapor<R, const N: usize>(
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
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
impl FiniteMoistureSolution {
    /// Complete a duration through the original equal-pressure solution/vapor
    /// activity bridge. Nonvolatile solute amounts stay fixed; the same finite
    /// solvent inventory evolves along both fine half steps. Local errors check
    /// solvent mass even for pure solvent whose activity remains exactly one.
    pub fn advance_adaptive_to(
        &self,
        moisture: &MoistureSlab,
        target_s: f64,
        faces: &[Option<DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        physical: SolutionMoistureSolve,
        time: MoistureTimeControls,
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
        let initial = (self.clone(), moisture.clone());
        let ((reservoir, field), report) = evolve(
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
            cx,
        )?;
        Ok((reservoir, field, report))
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
        time: MoistureTimeControls,
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
        let initial = (self.clone(), moisture.clone());
        let ((reservoir, field), report) = evolve(
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
            cx,
        )?;
        Ok((reservoir, field, report))
    }
}

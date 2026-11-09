//! Two finite solution baths under the original conservative time controller.
//!
//! Each bath's activity, solvent mass and optional actual kinetic phase stocks
//! enter the comparison independently. Opposing exchanges must not disappear
//! into a combined reservoir coordinate. All nested field trials, including
//! rejected or failed temporal probes, share one explicit work allowance.

use super::*;
use crate::moisture::{
    radial::MoistureRadial,
    reservoir::water_sum,
    solution::{
        FiniteMoistureSolution, MoistureSolutionPort, PairedSolutionMoistureSolve,
        PairedSolutionMoistureStep, SolutionMoistureError, SolutionMoistureFace, advance_pair,
    },
};

/// Whole-duration transport accuracy, independent phase scales and actual root
/// work. The original per-endpoint root/field limits remain active as well.
#[derive(Clone, Copy, Debug)]
pub struct PairedSolutionMoistureTimeControls {
    pub transport: MoistureTimeControls,
    /// First and Second retain their own optional kinetic component tolerances.
    /// `Some` requires that bath's original finite-rate phase owner.
    pub kinetic_phases: [Option<KineticPhaseAccuracy>; 2],
    /// Cumulative native field proposals across every nested bath root and
    /// temporal trial, in 1..=100_000. Charged before surface/field admission.
    pub maximum_field_evaluations: usize,
}

/// One bath's original-to-final account, using only its assigned signed faces.
#[derive(Clone, Debug)]
pub struct SolutionMoistureAdvanceAccount {
    pub water_activity: [f64; 2],
    /// Thermodynamic diagnostic, never an applied mechanical surface traction.
    pub osmotic_pressure_pa: [f64; 2],
    pub solvent_mass_kg: [f64; 2],
    pub water_to_solid_kg: f64,
    pub water_balance_residual_kg: f64,
}

/// Both finite baths and the original native field over a completed duration.
#[derive(Clone, Debug)]
pub struct PairedSolutionMoistureAdvanceReport<Report = MoistureStepReport, const FACES: usize = 2>
{
    /// No single activity is assigned to two independent solutions. Reservoir
    /// water endpoints here are the sums of the two original solvent stores.
    pub transport: MoistureAdvanceReport<PairedSolutionMoistureStep<Report, FACES>, FACES>,
    pub solutions: [SolutionMoistureAdvanceAccount; 2],
    /// Actual nested field proposals; `transport.evaluations` counts endpoints.
    pub field_evaluations: usize,
    /// Independent sum of represented bath and cell endpoint/start inventories.
    pub direct_water_balance_residual_kg: f64,
}

impl<R: WaterStep<N>, const N: usize> WaterStep<N> for PairedSolutionMoistureStep<R, N> {
    fn boundary_water(&self) -> [f64; N] {
        self.transport.boundary_water()
    }
}

impl PairedSolutionMoistureTimeControls {
    pub(super) fn bath(self, port: usize) -> SolutionMoistureTimeControls {
        SolutionMoistureTimeControls {
            transport: self.transport,
            kinetic_phases: self.kinetic_phases[port],
        }
    }
}

pub(super) fn bath_error(
    old: &[FiniteMoistureSolution; 2],
    coarse: &[FiniteMoistureSolution; 2],
    fine: &[FiniteMoistureSolution; 2],
    controls: PairedSolutionMoistureTimeControls,
    cx: &Cx<'_>,
) -> Result<f64, MoistureAdvanceError<SolutionMoistureError>> {
    let mut ratio = 0.0_f64;
    for port in 0..2 {
        checkpoint(cx)?;
        let before = &old[port];
        let a = &coarse[port];
        let b = &fine[port];
        let activity = [before, a, b].map(|bath| {
            bath.solution()
                .water_activity()
                .map_err(SolutionMoistureError::from)
                .map_err(MoistureAdvanceError::FinalInventory)
        });
        let [old_activity, coarse_activity, fine_activity] = activity;
        ratio = ratio.max(scaled_error(
            old_activity?,
            coarse_activity?,
            fine_activity?,
            controls.transport.activity_absolute_tolerance,
            controls.transport.relative_tolerance,
        )?);
        ratio = ratio.max(scaled_error(
            before.solvent_mass_kg(),
            a.solvent_mass_kg(),
            b.solvent_mass_kg(),
            controls.transport.reservoir_water_absolute_tolerance_kg,
            controls.transport.relative_tolerance,
        )?);
        ratio = ratio.max(controls.bath(port).phase_error(before, a, b, cx)?);
    }
    Ok(ratio)
}

pub(super) fn account(
    old: &FiniteMoistureSolution,
    next: &FiniteMoistureSolution,
    received: f64,
    tolerance: f64,
    cx: &Cx<'_>,
) -> Result<SolutionMoistureAdvanceAccount, MoistureAdvanceError<SolutionMoistureError>> {
    old.admit_exchanged_components(next, cx)
        .map_err(MoistureAdvanceError::FinalInventory)?;
    let solvent_mass_kg = [old.solvent_mass_kg(), next.solvent_mass_kg()];
    let residual = water_sum([solvent_mass_kg[1], -solvent_mass_kg[0], received]);
    if !residual.is_finite() || residual.abs() > tolerance {
        return Err(MoistureAdvanceError::WaterBalance {
            residual_kg: residual,
            tolerance_kg: tolerance,
        });
    }
    let [a, b] = [old, next].map(|bath| {
        bath.solution()
            .water_activity()
            .map_err(SolutionMoistureError::from)
            .map_err(MoistureAdvanceError::FinalInventory)
    });
    let [pa, pb] = [old, next].map(|bath| {
        bath.solution()
            .osmotic_pressure_pa()
            .map_err(SolutionMoistureError::from)
            .map_err(MoistureAdvanceError::FinalInventory)
    });
    Ok(SolutionMoistureAdvanceAccount {
        water_activity: [a?, b?],
        osmotic_pressure_pa: [pa?, pb?],
        solvent_mass_kg,
        water_to_solid_kg: received,
        water_balance_residual_kg: residual,
    })
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn advance<F: ReservoirField<N> + Clone, const N: usize>(
    original: [&FiniteMoistureSolution; 2],
    field: &F,
    target_s: f64,
    faces: &[Option<SolutionMoistureFace>; N],
    transport: F::Controls,
    physical: PairedSolutionMoistureSolve,
    controls: PairedSolutionMoistureTimeControls,
    field_error: impl Fn(&F, &F, &F) -> Result<f64, MoistureAdvanceError<SolutionMoistureError>>,
    cx: &Cx<'_>,
) -> Result<
    (
        [FiniteMoistureSolution; 2],
        F,
        PairedSolutionMoistureAdvanceReport<F::Report, N>,
    ),
    MoistureAdvanceError<SolutionMoistureError>,
>
where
    F::Report: WaterStep<N>,
{
    checkpoint(cx)?;
    if !(1..=100_000).contains(&controls.maximum_field_evaluations) {
        return Err(MoistureAdvanceError::Invalid(
            "paired solution durations need bounded cumulative native field work",
        ));
    }
    for (port, bath) in original.iter().enumerate() {
        controls.bath(port).admit(bath, cx)?;
    }
    let initial = ([original[0].clone(), original[1].clone()], field.clone());
    let mut field_evaluations = 0;
    let ((solutions, moisture), mut report) = evolve_with_error(
        &initial,
        target_s,
        controls.transport,
        |state, dt| {
            let (baths, field, step) = advance_pair(
                [&state.0[0], &state.0[1]],
                &state.1,
                dt,
                faces,
                transport,
                physical,
                cx,
                WaterStep::boundary_water,
                || {
                    if field_evaluations == controls.maximum_field_evaluations {
                        return Err(SolutionMoistureError::EvaluationBudget {
                            evaluations: field_evaluations,
                            maximum: controls.maximum_field_evaluations,
                        });
                    }
                    field_evaluations += 1;
                    Ok(())
                },
            )?
            .into_parts();
            Ok(((baths, field), step))
        },
        |state| field_snapshot::<F, N>(&state.1).map_err(SolutionMoistureError::from),
        |old, a, b| {
            Ok(bath_error(&old.0, &a.0, &b.0, controls, cx)?.max(field_error(&old.1, &a.1, &b.1)?))
        },
        cx,
    )?;
    let mut received = [WaterSum::default(); 2];
    for (side, face) in faces.iter().enumerate() {
        if let Some(face) = face {
            let port = usize::from(face.solution == MoistureSolutionPort::Second);
            received[port].add(report.boundary_water_kg[side]);
        }
    }
    let tolerance = controls.transport.water_balance_tolerance_kg;
    let accounts = [
        account(
            original[0],
            &solutions[0],
            received[0].total(),
            tolerance,
            cx,
        )?,
        account(
            original[1],
            &solutions[1],
            received[1].total(),
            tolerance,
            cx,
        )?,
    ];
    let old = field.inventory();
    let next = moisture.inventory();
    let bath_changes = [
        accounts[0].solvent_mass_kg[1],
        -accounts[0].solvent_mass_kg[0],
        accounts[1].solvent_mass_kg[1],
        -accounts[1].solvent_mass_kg[0],
    ];
    let balance = water_sum(bath_changes.into_iter().chain([report.water_change_kg]));
    let mut direct = WaterSum::default();
    for change in bath_changes {
        direct.add(change);
    }
    for (i, (&a, &b)) in old.moisture.iter().zip(next.moisture).enumerate() {
        if i % 128 == 0 {
            checkpoint(cx)?;
        }
        direct.add(b * next.dry_mass_at(i));
        direct.add(-a * old.dry_mass_at(i));
    }
    let direct = direct.total();
    for residual in [balance, direct] {
        if !residual.is_finite() || residual.abs() > tolerance {
            return Err(MoistureAdvanceError::WaterBalance {
                residual_kg: residual,
                tolerance_kg: tolerance,
            });
        }
    }
    let totals = std::array::from_fn(|i| {
        water_sum([
            accounts[0].solvent_mass_kg[i],
            accounts[1].solvent_mass_kg[i],
        ])
    });
    if totals.iter().any(|value| !value.is_finite()) {
        return Err(MoistureAdvanceError::Invalid(
            "combined solvent inventory is unrepresentable",
        ));
    }
    report.reservoir_water_kg = Some(totals);
    report.water_balance_residual_kg = balance;
    checkpoint(cx)?;
    Ok((
        solutions,
        moisture,
        PairedSolutionMoistureAdvanceReport {
            transport: report,
            solutions: accounts,
            field_evaluations,
            direct_water_balance_residual_kg: direct,
        },
    ))
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
impl FiniteMoistureSolution {
    /// Evolve the original slab and two finite baths to a complete target.
    /// Each bath's activity, water and requested phase stocks are resolved
    /// independently. Only accepted half-step transfers enter the water account.
    /// Full/half comparisons estimate local first-order error; no extrapolated
    /// state, second-order claim or global accuracy guarantee is constructed.
    /// A late source, work, balance or cancellation refusal returns no prefix.
    pub fn advance_paired_adaptive_to(
        &self,
        second: &Self,
        moisture: &MoistureSlab,
        target_s: f64,
        faces: &[Option<SolutionMoistureFace>; 2],
        transport: MoistureStepConfig,
        physical: PairedSolutionMoistureSolve,
        time: PairedSolutionMoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        ([Self; 2], MoistureSlab, PairedSolutionMoistureAdvanceReport),
        MoistureAdvanceError<SolutionMoistureError>,
    > {
        advance(
            [self, second],
            moisture,
            target_s,
            faces,
            transport,
            physical,
            time,
            |_, _, _| Ok(0.0),
            cx,
        )
    }

    /// The same bounded duration on all six original Cartesian faces. Source
    /// anisotropy, nonlinear diffusivity and physical surface resistances remain
    /// with the native field. Temperature and face assignments stay prescribed.
    pub fn advance_paired_cartesian_adaptive_to(
        &self,
        second: &Self,
        moisture: &MoistureCartesian,
        target_s: f64,
        faces: &[Option<SolutionMoistureFace>; 6],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: PairedSolutionMoistureSolve,
        time: PairedSolutionMoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            [Self; 2],
            MoistureCartesian,
            PairedSolutionMoistureAdvanceReport<MoistureCartesianReport, 6>,
        ),
        MoistureAdvanceError<SolutionMoistureError>,
    > {
        advance(
            [self, second],
            moisture,
            target_s,
            faces,
            (transport, nonlinear),
            physical,
            time,
            |_, _, _| Ok(0.0),
            cx,
        )
    }

    /// Cylindrical/spherical transport retains actual shell capacities and
    /// curved resistances. The existing radial L1 water-error norm supplements
    /// cell content and both baths' independent temporal coordinates.
    pub fn advance_paired_radial_adaptive_to(
        &self,
        second: &Self,
        moisture: &MoistureRadial,
        target_s: f64,
        faces: &[Option<SolutionMoistureFace>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: PairedSolutionMoistureSolve,
        time: PairedSolutionMoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            [Self; 2],
            MoistureRadial,
            PairedSolutionMoistureAdvanceReport,
        ),
        MoistureAdvanceError<SolutionMoistureError>,
    > {
        advance(
            [self, second],
            moisture,
            target_s,
            faces,
            (transport, nonlinear),
            physical,
            time,
            |old, a, b| radial::water_distance(old, a, b, time.transport, cx),
            cx,
        )
    }
}

//! Adaptive actual swelling geometry with two original finite solution baths.
//!
//! The same first-order endpoint owns transport, source mechanics and phase
//! kinetics in every full/half probe. Only conservative fine-path transfers
//! survive acceptance. Native displacement is a temporal observable; it does
//! not replace the separate nonlinear geometry or equilibrium tolerances.

use super::*;
use crate::moisture::{
    adaptive::solution_pair::{account, bath_error},
    reservoir::water_sum,
    solid::{
        MoistureLayeredSolid, MoistureLayeredSolidStep, MoistureSolidState,
        MoistureSpatialPairedSolutionControls, MoistureSpatialPairedSolutionReport,
    },
    solution::{
        FiniteMoistureSolution, MoistureSolutionPort, PairedSolutionMoistureStep,
        SolutionMoistureError,
    },
};
use std::sync::Arc;

/// Temporal controls for actual source-driven shape and both finite baths.
#[derive(Clone, Copy, Debug)]
pub struct MoistureSpatialPairedSolutionTimeControls {
    /// Local moisture/activity, dry-mass water, both actual solvent stores and
    /// optional independent dissolved/crystal stocks, with cumulative field work.
    pub solutions: PairedSolutionMoistureTimeControls,
    /// Absolute local error scale of each compatible nodal displacement [m].
    /// The shared relative tolerance multiplies displacement CHANGE from the
    /// last accepted state, never the reference coordinate or accumulated shape.
    pub displacement_absolute_tolerance_m: f64,
}

/// Complete-duration physical accounts and actual work, including rejected probes.
#[derive(Clone, Debug)]
pub struct MoistureSpatialPairedSolutionAdvanceReport {
    pub transport: LayeredMoistureAdvanceReport<
        PairedSolutionMoistureStep<MoistureSpatialPairedSolutionReport>,
    >,
    pub solutions: [SolutionMoistureAdvanceAccount; 2],
    /// All nested transport proposals; the summary counts coupled endpoints.
    pub field_evaluations: usize,
    /// All original native solves, including the initial displacement datum.
    pub native_solid_solves: usize,
    pub direct_water_balance_residual_kg: f64,
}

impl LayeredTransport for MoistureSpatialPairedSolutionReport {
    fn transport(&self) -> &MoistureLayeredReport {
        &self.spatial.transport
    }
}
impl<R: LayeredTransport> LayeredTransport for PairedSolutionMoistureStep<R> {
    fn transport(&self) -> &MoistureLayeredReport {
        self.transport.transport()
    }
}

fn displacement_error(
    old: &MoistureSolidState,
    coarse: &MoistureSolidState,
    fine: &MoistureSolidState,
    absolute: f64,
    relative: f64,
    cx: &Cx<'_>,
) -> Result<f64, MoistureAdvanceError<SolutionMoistureError>> {
    let [before, a, b] = [old, coarse, fine].map(|s| s.solution().displacement_m());
    if before.is_empty() || before.len() != a.len() || before.len() != b.len() {
        return Err(MoistureAdvanceError::Invalid(
            "temporal shape comparison requires the original native node layout",
        ));
    }
    let mut ratio = 0.0_f64;
    for (i, ((old, a), b)) in before.iter().zip(a).zip(b).enumerate() {
        if i % 128 == 0 {
            checkpoint(cx)?;
        }
        for axis in 0..3 {
            ratio = ratio.max(scaled_error(
                old[axis], a[axis], b[axis], absolute, relative,
            )?);
        }
    }
    Ok(ratio)
}

impl MoistureLayeredSolid {
    /// Evolve actual compatible swelling and two independently finite baths to
    /// a target time with the original bounded first-order step-doubling owner.
    ///
    /// Every full/half endpoint closes the original nested solvent roots and
    /// native geometric diffusion feedback. The same original solid binding,
    /// reference mesh, dry skeleton and swelling datum survive every probe and
    /// accepted interval. No independently expanded or extrapolated shape is
    /// fed back into transport. The returned native stress, wet inertia and
    /// prestressed modal pencil belong to the final accepted moisture state.
    ///
    /// In addition to per-layer moisture/activity and L1 physical-water error,
    /// each bath retains independent solvent/activity and optional actual phase
    /// tolerances. Compatible nodal displacement has its own dimensional scale.
    /// These are local asymptotic temporal estimates, not global error bounds or
    /// replacements for spatial, mechanical, root or geometry admissibility.
    ///
    /// `physical.maximum_native_solid_solves` bounds the ENTIRE duration,
    /// including one initial native state solve and every rejected bath/geometry
    /// trial. The time controls separately bound all transport proposals and
    /// coupled endpoints. Bath, layer and whole-system original-to-final water
    /// audits use fixed whole-window budgets, independent of accepted step count.
    /// Source, budget, geometry and cancellation refusals return no prefix.
    ///
    /// Temperature, exterior laws and vessel geometry remain prescribed. The
    /// original equal-pressure vapor bridge, separable prism, scalar current-
    /// frame diffusivity and small-strain mechanical scope remain unchanged.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_spatial_diffusion_paired_solutions_adaptive_to(
        &self,
        solutions: [&FiniteMoistureSolution; 2],
        moisture: &MoistureLayeredSlab,
        target_s: f64,
        faces: [Option<MoistureSolutionPort>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: MoistureSpatialPairedSolutionControls,
        accuracy: MoistureSpatialPairedSolutionTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            [FiniteMoistureSolution; 2],
            MoistureLayeredSolidStep<MoistureSpatialPairedSolutionAdvanceReport>,
        ),
        MoistureAdvanceError<SolutionMoistureError>,
    > {
        checkpoint(cx)?;
        let time = accuracy.solutions.transport;
        validate(moisture.time_s(), target_s, time)?;
        if !positive(accuracy.displacement_absolute_tolerance_m)
            || !positive(physical.maximum_geometry_residual_m)
            || !(1..=100_000).contains(&physical.maximum_native_solid_solves)
            || !(1..=100_000).contains(&accuracy.solutions.maximum_field_evaluations)
        {
            return Err(MoistureAdvanceError::Invalid(
                "adaptive spatial exchange needs positive displacement/geometry scales and bounded actual work",
            ));
        }
        for (port, bath) in solutions.iter().enumerate() {
            accuracy.solutions.bath(port).admit(bath, cx)?;
        }
        // Reserve before resolving the initial displacement comparison datum.
        // Arc makes controller snapshots share, rather than duplicate, native
        // sparse assemblies. Every proposed endpoint still runs actual mechanics.
        let mut native_solid_solves = 1;
        let initial_solid =
            self.resolve(moisture, cx)
                .map_err(|error| MoistureAdvanceError::Step {
                    time_s: moisture.time_s(),
                    step_s: 0.0,
                    evaluations: 0,
                    error: SolutionMoistureError::from(error),
                })?;
        let initial = (
            [solutions[0].clone(), solutions[1].clone()],
            moisture.clone(),
            Arc::new(initial_solid),
        );
        let mut field_evaluations = 0;
        let ((baths, field, solid), report) = evolve_fields_with_error(
            &initial,
            target_s,
            time,
            |state, dt| {
                let (candidate, solid) = self.advance_spatial_diffusion_paired_solutions_budgeted(
                    [&state.0[0], &state.0[1]],
                    &state.1,
                    dt,
                    faces,
                    transport,
                    nonlinear,
                    physical,
                    &mut native_solid_solves,
                    || {
                        if field_evaluations == accuracy.solutions.maximum_field_evaluations {
                            return Err(SolutionMoistureError::EvaluationBudget {
                                evaluations: field_evaluations,
                                maximum: accuracy.solutions.maximum_field_evaluations,
                            });
                        }
                        field_evaluations += 1;
                        Ok(())
                    },
                    cx,
                )?;
                let (baths, field, step) = candidate.into_parts();
                Ok(((baths, field, Arc::new(solid)), LayeredStep(step)))
            },
            |state| Ok(observe(&state.1)),
            |old, coarse, fine| {
                let water = water_distance(
                    old.1
                        .layers()
                        .iter()
                        .zip(coarse.1.layers())
                        .zip(fine.1.layers())
                        .map(|((a, b), c)| (a.slab().field(), b.slab().field(), c.slab().field())),
                    time,
                    cx,
                )?;
                Ok(water
                    .max(bath_error(
                        &old.0,
                        &coarse.0,
                        &fine.0,
                        accuracy.solutions,
                        cx,
                    )?)
                    .max(displacement_error(
                        &old.2,
                        &coarse.2,
                        &fine.2,
                        accuracy.displacement_absolute_tolerance_m,
                        time.relative_tolerance,
                        cx,
                    )?))
            },
            cx,
        )?;
        let tolerance = time.water_balance_tolerance_kg;
        let mut report = finish(report, tolerance, cx)?;
        let mut received = [WaterSum::default(); 2];
        for (side, assignment) in report
            .summary
            .last_substep
            .face_solutions
            .iter()
            .enumerate()
        {
            if let Some(port) = assignment {
                received[port.index()].add(report.summary.boundary_water_kg[side]);
            }
        }
        let accounts = [
            account(solutions[0], &baths[0], received[0].total(), tolerance, cx)?,
            account(solutions[1], &baths[1], received[1].total(), tolerance, cx)?,
        ];
        let bath_changes = [
            accounts[0].solvent_mass_kg[1],
            -accounts[0].solvent_mass_kg[0],
            accounts[1].solvent_mass_kg[1],
            -accounts[1].solvent_mass_kg[0],
        ];
        let balance = water_sum(
            bath_changes
                .into_iter()
                .chain([report.summary.water_change_kg]),
        );
        let mut direct = WaterSum::default();
        for change in bath_changes {
            direct.add(change);
        }
        for (old, next) in moisture.layers().iter().zip(field.layers()) {
            checkpoint(cx)?;
            let (old, next) = (old.slab().field(), next.slab().field());
            for (i, (&a, &b)) in old.moisture.iter().zip(next.moisture).enumerate() {
                if i % 128 == 0 {
                    checkpoint(cx)?;
                }
                direct.add(b * next.dry_mass_at(i));
                direct.add(-a * old.dry_mass_at(i));
            }
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
        let totals: [f64; 2] = std::array::from_fn(|i| {
            water_sum([
                accounts[0].solvent_mass_kg[i],
                accounts[1].solvent_mass_kg[i],
            ])
        });
        if totals.iter().any(|value| !value.is_finite()) {
            return Err(MoistureAdvanceError::Invalid(
                "combined spatial-exchange solvent inventory is unrepresentable",
            ));
        }
        report.summary.reservoir_water_kg = Some(totals);
        report.summary.water_balance_residual_kg = balance;
        checkpoint(cx)?;
        Ok((
            baths,
            MoistureLayeredSolidStep {
                moisture: field,
                transport: MoistureSpatialPairedSolutionAdvanceReport {
                    transport: report,
                    solutions: accounts,
                    field_evaluations,
                    native_solid_solves,
                    direct_water_balance_residual_kg: direct,
                },
                solid: Arc::unwrap_or_clone(solid),
            },
        ))
    }
}

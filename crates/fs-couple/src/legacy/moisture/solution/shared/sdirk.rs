//! Two-stage diffusion driven by a finite solution's evolving water activity.
//!
//! Each stage solves native specimens and the same original solvent inventory
//! together. Phase equilibrium is algebraic at each trial volume. Kinetic
//! crystal stocks need their own RK stage equations and are explicitly refused.

use super::*;
use crate::moisture::cartesian::MoistureCartesianSdirk2Report;
use crate::moisture::reservoir::{CartesianSdirkStage, SlabSdirkStage};
use crate::moisture::sdirk::{GAMMA, MoistureSdirk2Report};
use crate::moisture::surface::ResolvedMoistureFace;

/// Original-to-final finite-solution account with both native slab stages.
pub type SlabSolutionSdirk2Step = SolutionMoistureStep<MoistureSdirk2Report>;
/// One finite solution and native slab accepted after both coupled stages.
pub type SlabSolutionSdirk2Candidate =
    SolutionMoistureCandidate<MoistureSlab, MoistureSdirk2Report>;
/// Multiple native slabs sharing one second-order finite solvent account.
pub type SharedSlabSolutionSdirk2Candidate =
    SharedSolutionMoistureCandidate<MoistureSlab, MoistureSdirk2Report>;

/// Original-to-final solution account with both native spatial stage receipts.
pub type CartesianSolutionSdirk2Step = SolutionMoistureStep<MoistureCartesianSdirk2Report>;
/// One accepted solution and Cartesian field after both coupled stages.
pub type CartesianSolutionSdirk2Candidate =
    SolutionMoistureCandidate<MoistureCartesian, MoistureCartesianSdirk2Report>;
/// Multiple Cartesian fields sharing the same two-stage finite solution.
pub type SharedCartesianSolutionSdirk2Candidate =
    SharedSolutionMoistureCandidate<MoistureCartesian, MoistureCartesianSdirk2Report>;

impl FiniteMoistureSolution {
    /// Advance the native slab with second-order diffusion driven by this
    /// finite solution's actual activity. Uptake concentrates the solution;
    /// desorption dilutes it while nonvolatile component totals stay fixed.
    ///
    /// Both stages use the original source-owned solvent-volume root, sorption,
    /// surface exchange and slab diffusion, including Kirchhoff transport when
    /// selected. The endpoint includes each weighted transfer once; its first
    /// stage is diagnostic. The slab's ordinary integration mode is preserved.
    ///
    /// Neutral/electrolyte equilibrium remains algebraic at every trial volume.
    /// Finite-rate phase kinetics is refused until its actual stocks have
    /// matching RK stage equations. Accuracy assumes a smooth source/phase
    /// branch; equal prescribed temperature/pressure and negligible vapor
    /// storage retain the existing solution coupling's physical scope.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_sdirk2(
        &self,
        moisture: &MoistureSlab,
        dt_s: f64,
        faces: &[Option<DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SlabSolutionSdirk2Candidate, SolutionMoistureError> {
        let participant = MoistureParticipant {
            id: 0,
            slab: moisture,
            faces,
            transport,
        };
        let (reservoir, mut solids, shared) = self
            .advance_shared_sdirk2(
                &[participant],
                dt_s,
                moisture.moisture().len(),
                nonlinear,
                controls,
                cx,
            )?
            .into_parts();
        let solid = solids.pop().expect("one admitted solution slab");
        let step = SolutionMoistureStep {
            transport: solid.transport,
            water_activity: shared.water_activity,
            osmotic_pressure_pa: shared.osmotic_pressure_pa,
            solvent_mass_kg: shared.solvent_mass_kg,
            water_to_solid_kg: shared.water_to_solids_kg,
            water_balance_residual_kg: shared.water_balance_residual_kg,
            iterations: shared.iterations,
        };
        Ok(SolutionMoistureCandidate {
            reservoir,
            moisture: solid.moisture,
            step,
        })
    }

    /// Couple original slabs to one finite solution activity at each SDIRK
    /// stage. Each material, geometry, surface law and signed water transfer
    /// remains local; simultaneous uptake and drying share one solvent debit.
    ///
    /// Common nonlinear controls apply to every specimen while its original
    /// ordinary integration mode remains unchanged. Stable IDs, accepted clocks,
    /// temperature, source domains and aggregate cells retain the native
    /// solution owner's admission. Stage two starts from ORIGINAL solvent and
    /// solid storage and adds the actual first-stage face load exactly once.
    ///
    /// Source refusal, invalid stages, kinetics or cancellation publishes no
    /// partial bath or specimen. SDIRK stability does not imply an unconditional
    /// maximum principle; a shorter physical interval can be necessary.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_shared_sdirk2(
        &self,
        participants: &[MoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        nonlinear: MoistureNonlinearConfig,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedSlabSolutionSdirk2Candidate, SolutionMoistureError> {
        poll(cx)?;
        nonlinear.validate()?;
        admit_count(participants.len(), maximum_total_cells)?;
        if self.kinetic_solubility().is_some() {
            return Err(SolutionMoistureError::Invalid(
                "second-order solution coupling requires algebraic phase partitioning",
            ));
        }
        let mut ordered: Vec<_> = participants.iter().collect();
        ordered.sort_unstable_by_key(|body| body.id);
        let start = ordered[0].slab.time_s();
        let implicit_dt_s = GAMMA * dt_s;
        let first_time = start + implicit_dt_s;
        let end = start + dt_s;
        if !positive(dt_s)
            || !positive(implicit_dt_s)
            || !end.is_finite()
            || first_time <= start
            || first_time >= end
        {
            return Err(MoistureError::Invalid(
                "finite-solution slab SDIRK needs two representable forward stage times",
            )
            .into());
        }

        // Bound all copied state first; configure each specimen only once,
        // outside the scalar volume trials.
        let mut cells = 0usize;
        for body in &ordered {
            poll(cx)?;
            cells = cells
                .checked_add(body.slab.moisture().len())
                .filter(|count| *count <= maximum_total_cells)
                .ok_or(SolutionMoistureError::Invalid(
                    "shared-solution total cell budget exhausted",
                ))?;
        }
        let mut staged = Vec::with_capacity(ordered.len());
        for body in &ordered {
            poll(cx)?;
            let mut field = body.slab.clone();
            field.implicit_diffusivity = Some(nonlinear);
            if field.kirchhoff_diffusivity.is_some() {
                field.kirchhoff_diffusivity = Some(nonlinear);
            }
            staged.push(field);
        }
        let first_fields: Vec<_> = ordered
            .iter()
            .zip(&staged)
            .map(|(body, field)| ReservoirParticipant {
                id: body.id,
                field,
                faces: body.faces,
                transport: body.transport,
            })
            .collect();
        let (_, first_solids, first_report) = self
            .advance_fields(
                &first_fields,
                implicit_dt_s,
                maximum_total_cells,
                controls,
                cx,
            )?
            .into_parts();
        let temperature = self.solution().temperature_k();
        let first_activity = first_report.water_activity[1];
        let mut loads = Vec::with_capacity(ordered.len());
        for ((body, field), first) in ordered.iter().zip(&staged).zip(&first_solids) {
            poll(cx)?;
            debug_assert_eq!(body.id, first.id);
            let mut boundaries = std::array::from_fn(|_| ResolvedMoistureFace::default());
            for (side, law) in body.faces.iter().enumerate() {
                if let Some(law) = law {
                    boundaries[side] = law.resolve(temperature, first_activity)?;
                }
            }
            loads.push(field.sdirk_load(
                dt_s,
                &first.moisture,
                &first.transport,
                &boundaries,
                cx,
            )?);
        }
        // Preserve the native receipts but release diagnostic first-stage
        // fields before the second coupled search.
        let first_reports: Vec<_> = first_solids
            .into_iter()
            .map(|first| (first.id, first.transport))
            .collect();
        let stages: Vec<_> = staged
            .iter()
            .zip(&loads)
            .map(|(field, known)| SlabSdirkStage::new(field, known, implicit_dt_s))
            .collect();
        let last_fields: Vec<_> = ordered
            .iter()
            .zip(&stages)
            .map(|(body, field)| ReservoirParticipant {
                id: body.id,
                field,
                faces: body.faces,
                transport: body.transport,
            })
            .collect();

        // The first-stage solution is diagnostic, never a second solvent debit.
        let (reservoir, solids, mut step) = self
            .advance_fields(&last_fields, dt_s, maximum_total_cells, controls, cx)?
            .into_parts();
        step.iterations += first_report.iterations;
        let solids = solids
            .into_iter()
            .zip(first_reports)
            .zip(&ordered)
            .map(|((last, (first_id, first_transport)), original)| {
                debug_assert_eq!(last.id, first_id);
                debug_assert_eq!(last.id, original.id);
                let mut moisture = last.moisture.into_field();
                moisture.implicit_diffusivity = original.slab.implicit_diffusivity;
                moisture.kirchhoff_diffusivity = original.slab.kirchhoff_diffusivity;
                SharedMoistureSolid {
                    id: last.id,
                    moisture,
                    transport: MoistureSdirk2Report {
                        endpoint: last.transport,
                        first_stage: first_transport,
                    },
                }
            })
            .collect();
        poll(cx)?;
        Ok(SharedSolutionMoistureCandidate {
            reservoir,
            solids,
            step,
        })
    }

    /// Couple second-order Cartesian diffusion to finite solution activity.
    ///
    /// Uptake concentrates the actual solution and desorption dilutes it.
    /// Every stage refreshes dissolved activity, sorption, surface resistance
    /// and principal diffusivity from the existing owners. Optional neutral or
    /// electrolyte equilibrium retains its original total component stocks.
    ///
    /// The endpoint receipt contains the full RK-weighted water exchange;
    /// the first-stage receipt is diagnostic and must not be debited again.
    /// Second-order accuracy assumes a smooth source/phase branch. Prescribed
    /// temperature, equal pressure and negligible vapor storage retain the
    /// original solution coupling scope; no immersion permeability, ion
    /// transport, thermal or mechanical feedback is introduced.
    ///
    /// Kinetic phase evolution is refused because the existing endpoint-frozen
    /// phase owner is first order. Inadmissible stages and unsupported source
    /// trials refuse the entire proposal and may require a shorter interval.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_cartesian_sdirk2(
        &self,
        moisture: &MoistureCartesian,
        dt_s: f64,
        faces: &[Option<DryBasisSorption>; 6],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<CartesianSolutionSdirk2Candidate, SolutionMoistureError> {
        let participant = CartesianMoistureParticipant {
            id: 0,
            moisture,
            faces,
            transport,
            nonlinear,
        };
        let (reservoir, mut solids, shared) = self
            .advance_shared_cartesian_sdirk2(
                &[participant],
                dt_s,
                moisture.moisture().len(),
                controls,
                cx,
            )?
            .into_parts();
        let solid = solids.pop().expect("one admitted solution participant");
        let step = SolutionMoistureStep {
            transport: solid.transport,
            water_activity: shared.water_activity,
            osmotic_pressure_pa: shared.osmotic_pressure_pa,
            solvent_mass_kg: shared.solvent_mass_kg,
            water_to_solid_kg: shared.water_to_solids_kg,
            water_balance_residual_kg: shared.water_balance_residual_kg,
            iterations: shared.iterations,
        };
        Ok(SolutionMoistureCandidate {
            reservoir,
            moisture: solid.moisture,
            step,
        })
    }

    /// Advance heterogeneous Cartesian specimens against one finite solution.
    ///
    /// Simultaneous adsorption and desorption share one compensated water
    /// balance and a common activity at each stage. Stable IDs, source support,
    /// accepted clocks, temperatures and aggregate cell admission follow the
    /// original shared solution owner. Stage two starts from the ORIGINAL
    /// solvent and specimen storage and applies the physical first-stage flux
    /// load exactly once. Scalar iteration counts include both coupled solves.
    pub fn advance_shared_cartesian_sdirk2(
        &self,
        participants: &[CartesianMoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedCartesianSolutionSdirk2Candidate, SolutionMoistureError> {
        poll(cx)?;
        admit_count(participants.len(), maximum_total_cells)?;
        if self.kinetic_solubility().is_some() {
            return Err(SolutionMoistureError::Invalid(
                "second-order solution coupling requires algebraic phase partitioning",
            ));
        }
        let start = participants[0].moisture.time_s();
        let implicit_dt_s = GAMMA * dt_s;
        let first_time = start + implicit_dt_s;
        let end = start + dt_s;
        if !positive(dt_s)
            || !positive(implicit_dt_s)
            || !end.is_finite()
            || first_time <= start
            || first_time >= end
        {
            return Err(MoistureError::Invalid(
                "finite-solution SDIRK needs two representable forward stage times",
            )
            .into());
        }
        let (_, first_solids, first_report) = self
            .advance_shared_cartesian(
                participants,
                implicit_dt_s,
                maximum_total_cells,
                controls,
                cx,
            )?
            .into_parts();
        let mut ordered: Vec<_> = participants.iter().collect();
        ordered.sort_unstable_by_key(|body| body.id);
        let temperature = self.solution().temperature_k();
        let first_activity = first_report.water_activity[1];
        let mut loads = Vec::with_capacity(ordered.len());
        for (body, first) in ordered.iter().zip(&first_solids) {
            poll(cx)?;
            debug_assert_eq!(body.id, first.id);
            let mut boundaries = std::array::from_fn(|_| ResolvedMoistureFace::default());
            for (side, law) in body.faces.iter().enumerate() {
                if let Some(law) = law {
                    boundaries[side] = law.resolve(temperature, first_activity)?;
                }
            }
            loads.push(body.moisture.sdirk_load(
                dt_s,
                &first.moisture,
                &first.transport,
                &boundaries,
                cx,
            )?);
        }
        let fields: Vec<_> = ordered
            .iter()
            .zip(&loads)
            .map(|(body, known)| CartesianSdirkStage::new(body.moisture, known, implicit_dt_s))
            .collect();
        let staged: Vec<_> = ordered
            .iter()
            .zip(&fields)
            .map(|(body, field)| ReservoirParticipant {
                id: body.id,
                field,
                faces: body.faces,
                transport: (body.transport, body.nonlinear),
            })
            .collect();

        // Original storage is deliberate: the stage-one solution is diagnostic.
        let (reservoir, solids, mut step) = self
            .advance_fields(&staged, dt_s, maximum_total_cells, controls, cx)?
            .into_parts();
        step.iterations += first_report.iterations;
        let solids = solids
            .into_iter()
            .zip(first_solids)
            .map(|(last, first)| {
                debug_assert_eq!(last.id, first.id);
                SharedMoistureSolid {
                    id: last.id,
                    moisture: last.moisture.into_field(),
                    transport: MoistureCartesianSdirk2Report {
                        endpoint: last.transport,
                        first_stage: first.transport,
                    },
                }
            })
            .collect();
        poll(cx)?;
        Ok(SharedSolutionMoistureCandidate {
            reservoir,
            solids,
            step,
        })
    }
}

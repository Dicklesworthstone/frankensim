//! Two-stage diffusion driven by a finite solution's evolving water activity.
//!
//! Each stage solves native specimens and the same original solvent inventory
//! together. Phase equilibrium is algebraic at each trial volume. Kinetic
//! crystal stocks need their own RK stage equations and are explicitly refused.

use super::*;
use crate::moisture::cartesian::MoistureCartesianSdirk2Report;
use crate::moisture::reservoir::CartesianSdirkStage;
use crate::moisture::sdirk::GAMMA;
use crate::moisture::surface::ResolvedMoistureFace;

/// Original-to-final solution account with both native spatial stage receipts.
pub type CartesianSolutionSdirk2Step = SolutionMoistureStep<MoistureCartesianSdirk2Report>;
/// One accepted solution and Cartesian field after both coupled stages.
pub type CartesianSolutionSdirk2Candidate =
    SolutionMoistureCandidate<MoistureCartesian, MoistureCartesianSdirk2Report>;
/// Multiple Cartesian fields sharing the same two-stage finite solution.
pub type SharedCartesianSolutionSdirk2Candidate =
    SharedSolutionMoistureCandidate<MoistureCartesian, MoistureCartesianSdirk2Report>;

impl FiniteMoistureSolution {
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

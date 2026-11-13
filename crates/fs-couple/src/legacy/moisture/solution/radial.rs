//! Curved specimens exchange water with the original finite aqueous inventory.
//!
//! Actual unequal shell dry masses enter the existing shared solution root.
//! Curvature, source-integrated diffusivity and surface resistance remain owned
//! by the native radial solver. This retains the equal-pressure equilibrated
//! vapor path: solution activity controls RH; it is not liquid immersion.

use super::*;
use crate::moisture::{
    MoistureSdirk2Report, radial::MoistureRadial, reservoir::RadialSdirkStage, sdirk::GAMMA,
    surface::ResolvedMoistureFace,
};

/// One finite solution and one original curved-shell moisture field.
pub type RadialSolutionMoistureCandidate = SolutionMoistureCandidate<MoistureRadial>;
/// Independently sourced radial fields sharing one finite solvent inventory.
pub type SharedRadialSolutionMoistureCandidate = SharedSolutionMoistureCandidate<MoistureRadial>;
/// Complete original-to-final solvent account with both native radial stages.
pub type RadialSolutionSdirk2Step = SolutionMoistureStep<MoistureSdirk2Report>;
/// One solution and radial field after conservative second-order integration.
pub type RadialSolutionSdirk2Candidate =
    SolutionMoistureCandidate<MoistureRadial, MoistureSdirk2Report>;
/// One solution and all radial specimens after both coupled stages.
pub type SharedRadialSolutionSdirk2Candidate =
    SharedSolutionMoistureCandidate<MoistureRadial, MoistureSdirk2Report>;

impl FiniteMoistureSolution {
    /// Exchange solvent with a cylinder, annulus, sphere or spherical shell.
    ///
    /// Uptake concentrates the original solution; desorption dilutes it. The
    /// retained solution owner resolves ideal/nonideal activity and any admitted
    /// finite phase partition at every solvent-volume trial. Nonvolatile solute
    /// remains there. `None` seals a face; a zero-area center must be sealed.
    /// No surrogate slab, uniform-cell inventory or liquid permeability is used.
    /// Temperature is prescribed and must match the original material sources.
    /// Both candidates remain private until source and water balances pass.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_radial(
        &self,
        moisture: &MoistureRadial,
        dt_s: f64,
        faces: &[Option<DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<RadialSolutionMoistureCandidate, SolutionMoistureError> {
        self.advance_field(moisture, dt_s, faces, (transport, nonlinear), controls, cx)
    }

    /// Solve one common endpoint activity for all original radial specimens.
    ///
    /// Simultaneous uptake and desorption retain their actual signed shell-water
    /// changes and debit the common solvent once. Stable distinct IDs determine
    /// summation order. Clocks, temperatures, aggregate cells and source support
    /// are admitted by the existing shared finite-solution owner. Sealed fields
    /// still undergo internal diffusion and selected reservoir phase evolution.
    pub fn advance_shared_radial(
        &self,
        participants: &[RadialMoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedRadialSolutionMoistureCandidate, SolutionMoistureError> {
        poll(cx)?;
        shared::admit_count(participants.len(), maximum_total_cells)?;
        let fields: Vec<_> = participants
            .iter()
            .map(|body| ReservoirParticipant {
                id: body.id,
                field: body.moisture,
                faces: body.faces,
                transport: (body.transport, body.nonlinear),
            })
            .collect();
        self.advance_fields(&fields, dt_s, maximum_total_cells, controls, cx)
    }

    /// Couple second-order radial diffusion to finite solution activity.
    ///
    /// Each stage queries its own dissolved activity, sorption and surface law.
    /// The endpoint account already contains the complete RK-weighted exchange;
    /// the diagnostic first stage must not be debited a second time. Invalid
    /// stages refuse transactionally and may require a shorter interval.
    /// Second-order accuracy assumes a smooth constitutive and phase branch.
    /// Kinetic crystal evolution is refused: its existing endpoint-frozen
    /// first-order law does not provide the required RK phase stages.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_radial_sdirk2(
        &self,
        moisture: &MoistureRadial,
        dt_s: f64,
        faces: &[Option<DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<RadialSolutionSdirk2Candidate, SolutionMoistureError> {
        let participant = RadialMoistureParticipant {
            id: 0,
            moisture,
            faces,
            transport,
            nonlinear,
        };
        let (reservoir, mut solids, shared) = self
            .advance_shared_radial_sdirk2(
                &[participant],
                dt_s,
                moisture.moisture().len(),
                controls,
                cx,
            )?
            .into_parts();
        let solid = solids
            .pop()
            .expect("one admitted radial solution participant");
        let step = SolutionMoistureStep {
            transport: solid.transport,
            water_activity: shared.water_activity,
            osmotic_pressure_pa: shared.osmotic_pressure_pa,
            solvent_mass_kg: shared.solvent_mass_kg,
            water_to_solid_kg: shared.water_to_solids_kg,
            water_balance_residual_kg: shared.water_balance_residual_kg,
            iterations: shared.iterations,
        };
        poll(cx)?;
        Ok(SolutionMoistureCandidate {
            reservoir,
            moisture: solid.moisture,
            step,
        })
    }

    /// Advance all radial specimens against the same finite solution at each
    /// SDIRK stage. Both stage roots start from the ORIGINAL solvent inventory.
    /// Stage two retains the actual first-stage curved-face flux load exactly
    /// once, normalized by each field's original shell storage measure.
    pub fn advance_shared_radial_sdirk2(
        &self,
        participants: &[RadialMoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedRadialSolutionSdirk2Candidate, SolutionMoistureError> {
        poll(cx)?;
        shared::admit_count(participants.len(), maximum_total_cells)?;
        if self.kinetic_solubility().is_some() {
            return Err(SolutionMoistureError::Invalid(
                "second-order radial solution coupling requires algebraic phase partitioning",
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
                "finite-solution radial SDIRK needs two representable forward stage times",
            )
            .into());
        }
        let (_, first_solids, first_report) = self
            .advance_shared_radial(
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
            loads.push(
                body.moisture
                    .sdirk_load(dt_s, &first.moisture, &boundaries, cx)?,
            );
        }
        let fields: Vec<_> = ordered
            .iter()
            .zip(&loads)
            .map(|(body, known)| RadialSdirkStage::new(body.moisture, known, implicit_dt_s))
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
                    transport: MoistureSdirk2Report {
                        endpoint: last.transport,
                        first_stage: first.transport,
                    },
                }
            })
            .collect();
        poll(cx)?;
        Ok(SharedSolutionMoistureCandidate::from_accepted(
            reservoir, solids, step,
        ))
    }
}

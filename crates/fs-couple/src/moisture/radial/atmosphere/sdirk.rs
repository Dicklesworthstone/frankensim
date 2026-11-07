//! Two-stage radial diffusion with the actual finite humidity at each stage.
//!
//! Stage one resolves the shared gas and every native radial field together.
//! Its accepted RH supplies the first-stage face rates. Stage two adds their
//! `(1-gamma)*dt` load to the ORIGINAL solid storage and balances the ORIGINAL
//! gas inventory against the final weighted shell-water increment. The first
//! gas endpoint is diagnostic; it is never debited as a second initial store.

use super::{RadialMoistureParticipant, water_change};
use crate::moisture::{
    MoistureError, MoistureNonlinearConfig, MoistureSdirk2Report, MoistureStepConfig,
    atmosphere::{
        AtmosphereMoistureCandidate, AtmosphereMoistureError, AtmosphereMoistureSolve,
        AtmosphereMoistureStep, DryBasisSorption, FiniteMoistureAtmosphere,
        SharedAtmosphereMoistureCandidate, SharedMoistureSolid, admit_participants,
    },
    poll, positive,
    radial::MoistureRadial,
    reservoir::water_sum,
    sdirk::GAMMA,
    surface::ResolvedMoistureFace,
};
use fs_exec::Cx;

/// Original-to-final finite vapor accounting with both native radial stages.
pub type RadialAtmosphereSdirk2Step = AtmosphereMoistureStep<MoistureSdirk2Report>;
/// One gas inventory and one radial field accepted after both coupled stages.
pub type RadialAtmosphereSdirk2Candidate =
    AtmosphereMoistureCandidate<MoistureRadial, MoistureSdirk2Report>;
/// Independently sourced radial fields sharing one second-order vapor account.
pub type SharedRadialAtmosphereSdirk2Candidate =
    SharedAtmosphereMoistureCandidate<MoistureRadial, MoistureSdirk2Report>;

fn boundaries(
    laws: &[Option<DryBasisSorption>; 2],
    temperature: f64,
    rh: f64,
    cx: &Cx<'_>,
) -> Result<[ResolvedMoistureFace; 2], AtmosphereMoistureError> {
    let mut result = std::array::from_fn(|_| ResolvedMoistureFace::default());
    for (side, law) in laws.iter().enumerate() {
        poll(cx)?;
        if let Some(law) = law {
            result[side] = law.resolve(temperature, rh)?;
        }
    }
    Ok(result)
}

impl FiniteMoistureAtmosphere {
    /// Second-order radial diffusion coupled to one actual finite atmosphere.
    ///
    /// Sorption, finite surface exchange and nonlinear diffusivity are evaluated
    /// at each stage's own RH. The endpoint report already contains the weighted
    /// water transfers of both stages; its diagnostic first-stage transfer must
    /// not be added again. All original source, range and mass checks remain.
    /// An inadmissible stage refuses the whole interval without clipping or a
    /// fallback to backward Euler; a shorter interval may be required.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_radial_sdirk2(
        &self,
        moisture: &MoistureRadial,
        dt_s: f64,
        faces: &[Option<DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<RadialAtmosphereSdirk2Candidate, AtmosphereMoistureError> {
        let participant = RadialMoistureParticipant {
            id: 0,
            moisture,
            faces,
            transport,
            nonlinear,
        };
        let (atmosphere, mut solids, shared) = self
            .advance_shared_radial_sdirk2(
                &[participant],
                dt_s,
                moisture.moisture().len(),
                controls,
                cx,
            )?
            .into_parts();
        let solid = solids.pop().expect("one admitted radial specimen");
        let step = AtmosphereMoistureStep {
            sorption: solid.transport.endpoint.surface_sorption.clone(),
            transport: solid.transport,
            relative_humidity: shared.relative_humidity,
            vapor_mass_kg: shared.vapor_mass_kg,
            water_to_solid_kg: shared.water_to_solids_kg,
            water_balance_residual_kg: shared.water_balance_residual_kg,
            maximum_carrier_mass_residual_kg: shared.maximum_carrier_mass_residual_kg,
            gas_internal_energy_change_j: shared.gas_internal_energy_change_j,
            iterations: shared.iterations,
        };
        poll(cx)?;
        Ok(AtmosphereMoistureCandidate::from_accepted(
            atmosphere,
            solid.moisture,
            step,
        ))
    }

    /// Advance radial specimens with one shared finite humidity at each stage.
    ///
    /// Actual unequal shell masses weight the compensated water account. Each
    /// field retains its material, geometry and face laws; stable IDs, common
    /// clocks, temperature and total cell limits follow the ordinary shared
    /// radial owner. Both scalar solves start from the original atmosphere.
    /// Temperature, free gas volume and reference solid geometry stay fixed;
    /// this does not add sorption heat, condensation or a temporal certificate.
    pub fn advance_shared_radial_sdirk2(
        &self,
        participants: &[RadialMoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedRadialAtmosphereSdirk2Candidate, AtmosphereMoistureError> {
        poll(cx)?;
        admit_participants(participants.len(), maximum_total_cells)?;
        let start = participants[0].moisture.time_s();
        let implicit_dt = GAMMA * dt_s;
        let first_time = start + implicit_dt;
        let end = start + dt_s;
        if !positive(dt_s)
            || !positive(implicit_dt)
            || !end.is_finite()
            || first_time <= start
            || first_time >= end
        {
            return Err(MoistureError::Invalid(
                "finite-vapor radial SDIRK needs two representable forward stage times",
            )
            .into());
        }
        let (_, first_solids, first_report) = self
            .advance_shared_radial(participants, implicit_dt, maximum_total_cells, controls, cx)?
            .into_parts();
        let mut ordered: Vec<_> = participants.iter().collect();
        ordered.sort_unstable_by_key(|body| body.id);
        let temperature = self.gas().state().state().temperature;
        let first_rh = first_report.relative_humidity[1];
        let mut loads = Vec::with_capacity(ordered.len());
        for (body, first) in ordered.iter().zip(&first_solids) {
            poll(cx)?;
            debug_assert_eq!(body.id, first.id);
            let faces = boundaries(body.faces, temperature, first_rh, cx)?;
            loads.push(
                body.moisture
                    .sdirk_load(dt_s, &first.moisture, &faces, cx)?,
            );
        }
        let all_sealed = ordered
            .iter()
            .all(|body| body.faces.iter().all(Option::is_none));

        // The full-interval solve deliberately uses self, never the temporary
        // first-stage gas. The known face load is already in each native report.
        let (atmosphere, last_solids, mut step) =
            self.advance_exchange(start, all_sealed, controls, cx, |rh| {
                let mut solids = Vec::with_capacity(ordered.len());
                for (body, known) in ordered.iter().zip(&loads) {
                    poll(cx)?;
                    let faces = boundaries(body.faces, temperature, rh, cx)?;
                    let (moisture, transport) = body.moisture.advance_stage(
                        dt_s,
                        implicit_dt,
                        faces,
                        body.transport,
                        body.nonlinear,
                        Some(known),
                        cx,
                    )?;
                    solids.push(SharedMoistureSolid {
                        id: body.id,
                        moisture,
                        transport,
                    });
                }
                let change = water_sum(
                    ordered
                        .iter()
                        .zip(&solids)
                        .map(|(before, after)| water_change(before.moisture, &after.moisture)),
                );
                let end_time = solids[0].moisture.time_s();
                Ok((solids, change, end_time))
            })?;
        step.iterations += first_report.iterations;
        step.maximum_carrier_mass_residual_kg = step
            .maximum_carrier_mass_residual_kg
            .max(first_report.maximum_carrier_mass_residual_kg);
        let solids = last_solids
            .into_iter()
            .zip(first_solids)
            .map(|(last, first)| {
                debug_assert_eq!(last.id, first.id);
                SharedMoistureSolid {
                    id: last.id,
                    moisture: last.moisture,
                    transport: MoistureSdirk2Report {
                        endpoint: last.transport,
                        first_stage: first.transport,
                    },
                }
            })
            .collect();
        poll(cx)?;
        Ok(SharedAtmosphereMoistureCandidate::from_accepted(
            atmosphere, solids, step,
        ))
    }
}

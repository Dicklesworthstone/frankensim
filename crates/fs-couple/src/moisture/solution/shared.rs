//! One finite solution root for any admitted collection of native specimens.

use super::field::SolutionParticipant;
use super::*;
use crate::moisture::reservoir::water_sum;

mod sdirk;
pub use sdirk::{
    CartesianSolutionSdirk2Candidate, CartesianSolutionSdirk2Step,
    SharedCartesianSolutionSdirk2Candidate, SharedSlabSolutionSdirk2Candidate,
    SlabSolutionSdirk2Candidate, SlabSolutionSdirk2Step,
};

/// One common solution debit, with all signed specimen transfers retained in
/// their original native reports. Total solute amounts are unchanged, including
/// finite solid stocks when selected, even when opposing adsorption/desorption
/// leaves the net solvent inventory unchanged.
#[derive(Clone, Debug)]
pub struct SharedSolutionMoistureStep {
    pub interval_s: [f64; 2],
    pub water_activity: [f64; 2],
    /// Diagnostic only; this is not an applied mechanical pressure.
    pub osmotic_pressure_pa: [f64; 2],
    pub solvent_mass_kg: [f64; 2],
    pub solid_water_kg: [f64; 2],
    pub water_to_solids_kg: f64,
    pub water_balance_residual_kg: f64,
    /// Independent balance using original/final represented stored inventories.
    pub direct_water_balance_residual_kg: f64,
    /// Interior scalar trials; initial/bracket-endpoint trials are additional.
    pub iterations: usize,
}

/// One actual finite solution and every original field admitted atomically.
/// Specimens are in ascending stable-ID order, independent of enumeration.
#[derive(Clone, Debug)]
pub struct SharedSolutionMoistureCandidate<Field = MoistureSlab, Report = MoistureStepReport> {
    reservoir: FiniteMoistureSolution,
    solids: Vec<SharedMoistureSolid<Field, Report>>,
    step: SharedSolutionMoistureStep,
}
impl<Field, Report> SharedSolutionMoistureCandidate<Field, Report> {
    pub(super) fn from_accepted(
        reservoir: FiniteMoistureSolution,
        solids: Vec<SharedMoistureSolid<Field, Report>>,
        step: SharedSolutionMoistureStep,
    ) -> Self {
        Self {
            reservoir,
            solids,
            step,
        }
    }

    pub const fn reservoir(&self) -> &FiniteMoistureSolution {
        &self.reservoir
    }
    pub fn solids(&self) -> &[SharedMoistureSolid<Field, Report>] {
        &self.solids
    }
    pub const fn step(&self) -> &SharedSolutionMoistureStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        FiniteMoistureSolution,
        Vec<SharedMoistureSolid<Field, Report>>,
        SharedSolutionMoistureStep,
    ) {
        (self.reservoir, self.solids, self.step)
    }
}

/// One finite aqueous inventory and complete, independently sourced 3D fields.
pub type SharedCartesianSolutionMoistureCandidate =
    SharedSolutionMoistureCandidate<MoistureCartesian, MoistureCartesianReport>;

pub(super) fn admit_count(count: usize, maximum: usize) -> Result<(), SolutionMoistureError> {
    if count == 0
        || maximum == 0
        || maximum > crate::moisture::MAX_MOISTURE_CELLS
        || count > maximum
    {
        return Err(SolutionMoistureError::Invalid(
            "a nonempty collection of solution-coupled specimens within the aggregate cell budget is required",
        ));
    }
    Ok(())
}

impl FiniteMoistureSolution {
    /// Simultaneously exchange water with independently sourced slabs through
    /// one common endpoint solution activity. Every trial restarts all original
    /// fields; a specimen's uptake can be supplied by another's desorption.
    ///
    /// Stable IDs must be unique, accepted clocks and temperatures must agree,
    /// and the total cell budget includes sealed specimens whose internal
    /// diffusion still advances. The original scalar and native solver work
    /// bounds apply to the entire coupled proposal and each trial respectively.
    /// All participants and the one reservoir publish together or not at all.
    pub fn advance_shared(
        &self,
        participants: &[MoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedSolutionMoistureCandidate, SolutionMoistureError> {
        poll(cx)?;
        admit_count(participants.len(), maximum_total_cells)?;
        let fields: Vec<_> = participants
            .iter()
            .map(|body| ReservoirParticipant {
                id: body.id,
                field: body.slab,
                faces: body.faces,
                transport: body.transport,
            })
            .collect();
        self.advance_fields(&fields, dt_s, maximum_total_cells, controls, cx)
    }

    /// Simultaneous six-face exchange for genuine Cartesian fields. Different
    /// geometries, source diffusivities and local surface laws share one solvent
    /// inventory and one water-activity root, without sequential reservoir debits.
    pub fn advance_shared_cartesian(
        &self,
        participants: &[CartesianMoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedCartesianSolutionMoistureCandidate, SolutionMoistureError> {
        poll(cx)?;
        admit_count(participants.len(), maximum_total_cells)?;
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

    pub(in crate::moisture) fn advance_fields<F: ReservoirField<N>, const N: usize>(
        &self,
        participants: &[ReservoirParticipant<'_, F, N>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedSolutionMoistureCandidate<F, F::Report>, SolutionMoistureError> {
        self.advance_participants(participants, dt_s, maximum_total_cells, controls, cx)
    }

    pub(super) fn advance_participants<P: SolutionParticipant>(
        &self,
        participants: &[P],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedSolutionMoistureCandidate<P::Field, P::Report>, SolutionMoistureError> {
        poll(cx)?;
        admit_count(participants.len(), maximum_total_cells)?;
        let mut ordered: Vec<_> = participants.iter().collect();
        ordered.sort_unstable_by_key(|body| body.id());
        if ordered.windows(2).any(|pair| pair[0].id() == pair[1].id()) {
            return Err(SolutionMoistureError::Invalid(
                "shared-solution specimen IDs must be unique",
            ));
        }
        let time_s = ordered[0].time_s();
        let temperature = self.solution().temperature_k();
        let mut cells = 0usize;
        for body in &ordered {
            poll(cx)?;
            cells = cells
                .checked_add(body.cells())
                .filter(|count| *count <= maximum_total_cells)
                .ok_or(SolutionMoistureError::Invalid(
                    "shared-solution total cell budget exhausted",
                ))?;
            if body.time_s().to_bits() != time_s.to_bits()
                || body.temperature_kelvin().to_bits() != temperature.to_bits()
            {
                return Err(SolutionMoistureError::Invalid(
                    "shared solution specimens need the same accepted clock and solution temperature",
                ));
            }
            body.admit(cx)?;
        }
        let (reservoir, solids, exchange) = self.advance_exchange(
            dt_s,
            ordered.iter().all(|body| body.sealed()),
            controls,
            cx,
            |activity| {
                let mut solids = Vec::with_capacity(ordered.len());
                for body in &ordered {
                    poll(cx)?;
                    let (moisture, transport) = body.advance(dt_s, temperature, activity, cx)?;
                    solids.push(SharedMoistureSolid {
                        id: body.id(),
                        moisture,
                        transport,
                    });
                }
                // Each original dry inventory appears once, in stable ID order.
                let water = water_sum(
                    ordered
                        .iter()
                        .zip(&solids)
                        .map(|(body, next)| body.water_change(&next.moisture)),
                );
                Ok((solids, water))
            },
        )?;
        let solvent_mass_kg = exchange.solvent_mass_kg;
        let residual = exchange.water_balance_residual_kg;
        let mut originals = Vec::with_capacity(ordered.len());
        let mut endpoints = Vec::with_capacity(ordered.len());
        for (body, next) in ordered.iter().zip(&solids) {
            poll(cx)?;
            originals.push(body.water_mass_kg(None));
            endpoints.push(body.water_mass_kg(Some(&next.moisture)));
        }
        let solid_water_kg = [
            water_sum(originals.iter().copied()),
            water_sum(endpoints.iter().copied()),
        ];
        let direct_residual = water_sum(
            [solvent_mass_kg[1], -solvent_mass_kg[0]].into_iter().chain(
                endpoints
                    .iter()
                    .zip(&originals)
                    .flat_map(|(&new, &old)| [new, -old]),
            ),
        );
        if solvent_mass_kg
            .iter()
            .chain(&solid_water_kg)
            .any(|m| !m.is_finite())
        {
            return Err(SolutionMoistureError::Invalid(
                "shared solution original/final water inventory is unrepresentable",
            ));
        }
        for residual_kg in [residual, direct_residual] {
            if !residual_kg.is_finite() || residual_kg.abs() > controls.mass_tolerance_kg {
                return Err(SolutionMoistureError::MassBalance {
                    residual_kg,
                    tolerance_kg: controls.mass_tolerance_kg,
                });
            }
        }
        let step = SharedSolutionMoistureStep {
            interval_s: [time_s, time_s + dt_s],
            water_activity: exchange.water_activity,
            osmotic_pressure_pa: exchange.osmotic_pressure_pa,
            solvent_mass_kg,
            solid_water_kg,
            water_to_solids_kg: exchange.water_to_solid_kg,
            water_balance_residual_kg: residual,
            direct_water_balance_residual_kg: direct_residual,
            iterations: exchange.iterations,
        };
        poll(cx)?;
        Ok(SharedSolutionMoistureCandidate {
            reservoir,
            solids,
            step,
        })
    }
}

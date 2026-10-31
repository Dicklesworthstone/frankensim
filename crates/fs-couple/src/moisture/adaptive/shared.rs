//! Thin bindings from the original simultaneous reservoirs to temporal control.

use super::*;
use crate::moisture::{
    atmosphere::{
        AtmosphereMoistureError, AtmosphereMoistureSolve, CartesianMoistureParticipant,
        FiniteMoistureAtmosphere, MoistureParticipant, SharedAtmosphereMoistureStep,
        SharedMoistureSolid,
    },
    reservoir::ReservoirParticipant,
    solution::{
        FiniteMoistureSolution, SharedSolutionMoistureStep, SolutionMoistureError,
        SolutionMoistureSolve,
    },
};

/// The complete reservoir and native specimen endpoints admitted together.
/// Each solid retains its last accepted native transport/source receipt;
/// `report` contains complete-window accounts in the same stable-ID order.
#[derive(Clone, Debug)]
pub struct SharedMoistureAdvanceCandidate<Reservoir, Field, Transport, Step, const FACES: usize> {
    reservoir: Reservoir,
    solids: Vec<SharedMoistureSolid<Field, Transport>>,
    report: SharedMoistureAdvanceReport<Step, FACES>,
}
impl<B, F, R, T, const N: usize> SharedMoistureAdvanceCandidate<B, F, R, T, N> {
    pub const fn reservoir(&self) -> &B {
        &self.reservoir
    }
    pub fn solids(&self) -> &[SharedMoistureSolid<F, R>] {
        &self.solids
    }
    pub const fn report(&self) -> &SharedMoistureAdvanceReport<T, N> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        B,
        Vec<SharedMoistureSolid<F, R>>,
        SharedMoistureAdvanceReport<T, N>,
    ) {
        (self.reservoir, self.solids, self.report)
    }
}

pub type SharedAtmosphereMoistureAdvanceCandidate = SharedMoistureAdvanceCandidate<
    FiniteMoistureAtmosphere,
    MoistureSlab,
    MoistureStepReport,
    SharedAtmosphereMoistureStep,
    2,
>;
pub type SharedCartesianAtmosphereMoistureAdvanceCandidate = SharedMoistureAdvanceCandidate<
    FiniteMoistureAtmosphere,
    MoistureCartesian,
    MoistureCartesianReport,
    SharedAtmosphereMoistureStep,
    6,
>;
pub type SharedSolutionMoistureAdvanceCandidate = SharedMoistureAdvanceCandidate<
    FiniteMoistureSolution,
    MoistureSlab,
    MoistureStepReport,
    SharedSolutionMoistureStep,
    2,
>;
pub type SharedCartesianSolutionMoistureAdvanceCandidate = SharedMoistureAdvanceCandidate<
    FiniteMoistureSolution,
    MoistureCartesian,
    MoistureCartesianReport,
    SharedSolutionMoistureStep,
    6,
>;

struct SharedStep<Step, Transport> {
    reservoir: Step,
    transports: Vec<(u64, Transport)>,
}
impl<T, R: WaterStep<N>, const N: usize> ReportTransfers<N> for SharedStep<T, R> {
    fn transfers(&self) -> Specimens<(u64, [f64; N])> {
        Specimens::Many(
            self.transports
                .iter()
                .map(|(id, step)| (*id, step.boundary_water()))
                .collect(),
        )
    }
}

fn admit_count<E>(
    count: usize,
    maximum: usize,
    cx: &Cx<'_>,
) -> Result<(), MoistureAdvanceError<E>> {
    checkpoint(cx)?;
    if count == 0 || count > maximum || maximum > MAX_MOISTURE_CELLS {
        return Err(MoistureAdvanceError::Invalid(
            "a nonempty shared moisture collection within the aggregate cell bound is required",
        ));
    }
    Ok(())
}
fn slabs<'a, E>(
    participants: &[MoistureParticipant<'a>],
    maximum: usize,
    cx: &Cx<'_>,
) -> Result<Vec<ReservoirParticipant<'a, MoistureSlab, 2>>, MoistureAdvanceError<E>> {
    admit_count(participants.len(), maximum, cx)?;
    Ok(participants
        .iter()
        .map(|p| ReservoirParticipant {
            id: p.id,
            field: p.slab,
            faces: p.faces,
            transport: p.transport,
        })
        .collect())
}
fn cartesians<'a, E>(
    participants: &[CartesianMoistureParticipant<'a>],
    maximum: usize,
    cx: &Cx<'_>,
) -> Result<Vec<ReservoirParticipant<'a, MoistureCartesian, 6>>, MoistureAdvanceError<E>> {
    admit_count(participants.len(), maximum, cx)?;
    Ok(participants
        .iter()
        .map(|p| ReservoirParticipant {
            id: p.id,
            field: p.moisture,
            faces: p.faces,
            transport: (p.transport, p.nonlinear),
        })
        .collect())
}
fn observe_vapor(air: &FiniteMoistureAtmosphere) -> Result<[f64; 2], AtmosphereMoistureError> {
    Ok([air.vapor_mass_kg()?, air.relative_humidity()?])
}
fn observe_solution(bath: &FiniteMoistureSolution) -> Result<[f64; 2], SolutionMoistureError> {
    Ok([bath.solvent_mass_kg(), bath.solution().water_activity()?])
}
fn observe_fields<'a, B, F: ReservoirField<N>, E, const N: usize>(
    state: &'a (B, Vec<(u64, F)>),
    observe: fn(&B) -> Result<[f64; 2], E>,
) -> Result<Snapshot<'a>, E> {
    let [water, activity] = observe(&state.0)?;
    Ok(Snapshot {
        fields: Specimens::Many(
            state
                .1
                .iter()
                .map(|(id, field)| (*id, field.inventory()))
                .collect(),
        ),
        reservoir_water_kg: Some(water),
        activity: Some(activity),
    })
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn evolve_shared<
    B: Clone,
    F: ReservoirField<N> + Clone,
    T,
    E: Retryable + From<MoistureError>,
    const N: usize,
>(
    original: &B,
    mut participants: Vec<ReservoirParticipant<'_, F, N>>,
    target_s: f64,
    maximum_cells: usize,
    temperature: f64,
    time: MoistureTimeControls,
    observe: fn(&B) -> Result<[f64; 2], E>,
    mut propose: impl FnMut(
        &B,
        &[ReservoirParticipant<'_, F, N>],
        f64,
    ) -> Result<(B, Vec<SharedMoistureSolid<F, F::Report>>, T), E>,
    cx: &Cx<'_>,
) -> Result<SharedMoistureAdvanceCandidate<B, F, F::Report, T, N>, MoistureAdvanceError<E>>
where
    F::Report: WaterStep<N>,
{
    admit_count(participants.len(), maximum_cells, cx)?;
    participants.sort_unstable_by_key(|p| p.id);
    if participants.windows(2).any(|pair| pair[0].id == pair[1].id) {
        return Err(MoistureAdvanceError::Invalid(
            "shared adaptive specimen IDs must be unique",
        ));
    }
    let start = participants[0].field.inventory().time_s;
    validate(start, target_s, time)?;
    let mut cells = 0usize;
    for p in &participants {
        checkpoint(cx)?;
        let field = p.field.inventory();
        cells = cells
            .checked_add(field.moisture.len())
            .filter(|&n| n <= maximum_cells)
            .ok_or(MoistureAdvanceError::Invalid(
                "shared adaptive total cell budget exhausted",
            ))?;
        if field.time_s.to_bits() != start.to_bits()
            || field.reference.query().axes()["temperature"].to_bits() != temperature.to_bits()
        {
            return Err(MoistureAdvanceError::Invalid(
                "shared adaptive specimens need one accepted clock and reservoir temperature",
            ));
        }
    }
    let initial: (B, Vec<(u64, F)>) = (
        original.clone(),
        participants
            .iter()
            .map(|p| (p.id, (*p.field).clone()))
            .collect(),
    );
    let (state, report) = evolve_fields(
        &initial,
        target_s,
        time,
        |state, dt| {
            poll(cx)?;
            if state.1.len() != participants.len() {
                return Err(E::from(MoistureError::Invalid(
                    "shared native endpoint lost a specimen",
                )));
            }
            let mut fields = Vec::with_capacity(participants.len());
            for ((id, field), template) in state.1.iter().zip(&participants) {
                poll(cx)?;
                if *id != template.id {
                    return Err(E::from(MoistureError::Invalid(
                        "shared native endpoint changed specimen IDs",
                    )));
                }
                fields.push(ReservoirParticipant {
                    id: *id,
                    field,
                    faces: template.faces,
                    transport: template.transport,
                });
            }
            let (reservoir, solids, step) = propose(&state.0, &fields, dt)?;
            let mut next = Vec::with_capacity(solids.len());
            let mut transports = Vec::with_capacity(solids.len());
            for solid in solids {
                next.push((solid.id, solid.moisture));
                transports.push((solid.id, solid.transport));
            }
            Ok((
                (reservoir, next),
                SharedStep {
                    reservoir: step,
                    transports,
                },
            ))
        },
        |state| observe_fields::<B, F, E, N>(state, observe),
        cx,
    )?;
    let mut transports = Vec::new();
    let summary = report.summary.map_step(|step| {
        transports = step.transports;
        step.reservoir
    });
    if state.1.len() != transports.len()
        || state
            .1
            .iter()
            .zip(&transports)
            .any(|(field, step)| field.0 != step.0)
    {
        return Err(MoistureAdvanceError::Invalid(
            "final shared sources do not match the original specimen IDs",
        ));
    }
    let solids = state
        .1
        .into_iter()
        .zip(transports)
        .map(|((id, moisture), (_, transport))| SharedMoistureSolid {
            id,
            moisture,
            transport,
        })
        .collect();
    checkpoint(cx)?;
    Ok(SharedMoistureAdvanceCandidate {
        reservoir: state.0,
        solids,
        report: SharedMoistureAdvanceReport {
            summary,
            specimens: report.specimens,
        },
    })
}

#[allow(clippy::too_many_arguments)]
impl FiniteMoistureAtmosphere {
    /// Complete a duration for all original slabs and ONE finite vapor store.
    /// Every cell, RH and vapor inventory receives independent temporal error
    /// control. Stable-ID face accounts integrate only accepted fine steps;
    /// source, carrier, work or late conservation failure publishes nothing.
    pub fn advance_shared_adaptive_to(
        &self,
        participants: &[MoistureParticipant<'_>],
        target_s: f64,
        maximum_total_cells: usize,
        physical: AtmosphereMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        SharedAtmosphereMoistureAdvanceCandidate,
        MoistureAdvanceError<AtmosphereMoistureError>,
    > {
        let fields = slabs(participants, maximum_total_cells, cx)?;
        let mut candidate = evolve_shared(
            self,
            fields,
            target_s,
            maximum_total_cells,
            self.gas().state().state().temperature,
            time,
            observe_vapor,
            |air, fields, dt| {
                air.advance_fields(fields, dt, maximum_total_cells, physical, cx)
                    .map(|c| c.into_parts())
            },
            cx,
        )?;
        finite::check_vapor(
            self,
            &candidate.reservoir,
            candidate.report.summary.water_change_kg,
            time.water_balance_tolerance_kg,
            &mut candidate.report.summary,
            cx,
        )?;
        Ok(candidate)
    }

    /// The same complete-window vapor coupling for complete Cartesian fields.
    /// Each six-face ledger belongs to its own geometry and original dry mass.
    pub fn advance_shared_cartesian_adaptive_to(
        &self,
        participants: &[CartesianMoistureParticipant<'_>],
        target_s: f64,
        maximum_total_cells: usize,
        physical: AtmosphereMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        SharedCartesianAtmosphereMoistureAdvanceCandidate,
        MoistureAdvanceError<AtmosphereMoistureError>,
    > {
        let fields = cartesians(participants, maximum_total_cells, cx)?;
        let mut candidate = evolve_shared(
            self,
            fields,
            target_s,
            maximum_total_cells,
            self.gas().state().state().temperature,
            time,
            observe_vapor,
            |air, fields, dt| {
                air.advance_fields(fields, dt, maximum_total_cells, physical, cx)
                    .map(|c| c.into_parts())
            },
            cx,
        )?;
        finite::check_vapor(
            self,
            &candidate.reservoir,
            candidate.report.summary.water_change_kg,
            time.water_balance_tolerance_kg,
            &mut candidate.report.summary,
            cx,
        )?;
        Ok(candidate)
    }
}

#[allow(clippy::too_many_arguments)]
impl FiniteMoistureSolution {
    /// Complete simultaneous uptake/desorption through the original finite
    /// solution activity bridge, keeping every nonvolatile solute amount fixed.
    /// A zero net bath change never suppresses opposing specimen error estimates.
    pub fn advance_shared_adaptive_to(
        &self,
        participants: &[MoistureParticipant<'_>],
        target_s: f64,
        maximum_total_cells: usize,
        physical: SolutionMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<SharedSolutionMoistureAdvanceCandidate, MoistureAdvanceError<SolutionMoistureError>>
    {
        let fields = slabs(participants, maximum_total_cells, cx)?;
        evolve_shared(
            self,
            fields,
            target_s,
            maximum_total_cells,
            self.solution().temperature_k(),
            time,
            observe_solution,
            |bath, fields, dt| {
                bath.advance_fields(fields, dt, maximum_total_cells, physical, cx)
                    .map(|c| c.into_parts())
            },
            cx,
        )
    }

    /// Complete simultaneous Cartesian diffusion and finite-solution evolution.
    /// All original source, scalar-root and per-field solver bounds remain active
    /// within each counted probe. Final water tolerance never scales with work.
    pub fn advance_shared_cartesian_adaptive_to(
        &self,
        participants: &[CartesianMoistureParticipant<'_>],
        target_s: f64,
        maximum_total_cells: usize,
        physical: SolutionMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        SharedCartesianSolutionMoistureAdvanceCandidate,
        MoistureAdvanceError<SolutionMoistureError>,
    > {
        let fields = cartesians(participants, maximum_total_cells, cx)?;
        evolve_shared(
            self,
            fields,
            target_s,
            maximum_total_cells,
            self.solution().temperature_k(),
            time,
            observe_solution,
            |bath, fields, dt| {
                bath.advance_fields(fields, dt, maximum_total_cells, physical, cx)
                    .map(|c| c.into_parts())
            },
            cx,
        )
    }
}

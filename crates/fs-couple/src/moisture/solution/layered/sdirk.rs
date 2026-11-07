//! Native two-stage heterogeneous diffusion against one actual finite solution.
//!
//! First-stage bodies and bath solve together. Their converged physical face
//! rates are captured once. The final root starts from the ORIGINAL bath and
//! material storage, retaining those weighted rates while its own activity
//! drives the final implicit faces. Kinetic crystal histories need their own
//! stage equations and are refused; algebraic phase partition remains native.

use super::*;
use crate::moisture::layered::{
    MoistureLayeredSdirk2Report,
    kernel::{self, StageLoad},
};

/// Whole-interval finite-solution account plus both original material stages.
pub type LayeredSolutionSdirk2Step = SolutionMoistureStep<MoistureLayeredSdirk2Report>;
pub type LayeredSolutionSdirk2Candidate =
    SolutionMoistureCandidate<MoistureLayeredSlab, MoistureLayeredSdirk2Report>;
pub type SharedLayeredSolutionSdirk2Candidate =
    SharedSolutionMoistureCandidate<MoistureLayeredSlab, MoistureLayeredSdirk2Report>;
pub type LayeredRadialSolutionSdirk2Candidate =
    SolutionMoistureCandidate<MoistureLayeredRadial, MoistureLayeredSdirk2Report>;
pub type SharedLayeredRadialSolutionSdirk2Candidate =
    SharedSolutionMoistureCandidate<MoistureLayeredRadial, MoistureLayeredSdirk2Report>;

/// Only the native geometry-specific stage calls differ. Source admission and
/// inventories continue to come from the original finite-solution participant.
trait LayeredSdirkParticipant: SolutionParticipant<Report = MoistureLayeredReport> {
    fn stage_load(
        &self,
        dt_s: f64,
        first: &Self::Field,
        activity: f64,
        cx: &Cx<'_>,
    ) -> Result<StageLoad, SolutionMoistureError>;
    fn finish_stage(
        &self,
        dt_s: f64,
        known: &StageLoad,
        activity: f64,
        cx: &Cx<'_>,
    ) -> Result<(Self::Field, MoistureLayeredReport), SolutionMoistureError>;
}

impl LayeredSdirkParticipant for LayeredMoistureParticipant<'_> {
    fn stage_load(
        &self,
        dt_s: f64,
        first: &Self::Field,
        activity: f64,
        cx: &Cx<'_>,
    ) -> Result<StageLoad, SolutionMoistureError> {
        let faces = boundary(
            exposure(self.exposed, self.moisture.exterior_transfer()),
            activity,
        );
        Ok(self.moisture.sdirk_load(dt_s, first, &faces, cx)?)
    }
    fn finish_stage(
        &self,
        dt_s: f64,
        known: &StageLoad,
        activity: f64,
        cx: &Cx<'_>,
    ) -> Result<(Self::Field, MoistureLayeredReport), SolutionMoistureError> {
        let faces = boundary(
            exposure(self.exposed, self.moisture.exterior_transfer()),
            activity,
        );
        Ok(self.moisture.advance_sdirk_stage(
            dt_s,
            &faces,
            known,
            self.transport,
            self.nonlinear,
            cx,
        )?)
    }
}

impl LayeredSdirkParticipant for LayeredRadialMoistureParticipant<'_> {
    fn stage_load(
        &self,
        dt_s: f64,
        first: &Self::Field,
        activity: f64,
        cx: &Cx<'_>,
    ) -> Result<StageLoad, SolutionMoistureError> {
        let faces = boundary(
            exposure(self.exposed, self.moisture.exterior_transfer()),
            activity,
        );
        Ok(self.moisture.sdirk_load(dt_s, first, &faces, cx)?)
    }
    fn finish_stage(
        &self,
        dt_s: f64,
        known: &StageLoad,
        activity: f64,
        cx: &Cx<'_>,
    ) -> Result<(Self::Field, MoistureLayeredReport), SolutionMoistureError> {
        let faces = boundary(
            exposure(self.exposed, self.moisture.exterior_transfer()),
            activity,
        );
        Ok(self.moisture.advance_sdirk_stage(
            dt_s,
            &faces,
            known,
            self.transport,
            self.nonlinear,
            cx,
        )?)
    }
}

struct FinalLayeredStage<'a, P> {
    original: &'a P,
    known: &'a StageLoad,
}
impl<P: LayeredSdirkParticipant> SolutionParticipant for FinalLayeredStage<'_, P> {
    type Field = P::Field;
    type Report = MoistureLayeredReport;
    fn id(&self) -> u64 {
        self.original.id()
    }
    fn time_s(&self) -> f64 {
        self.original.time_s()
    }
    fn temperature_kelvin(&self) -> f64 {
        self.original.temperature_kelvin()
    }
    fn cells(&self) -> usize {
        self.original.cells()
    }
    fn sealed(&self) -> bool {
        self.original.sealed()
    }
    fn admit(&self, cx: &Cx<'_>) -> Result<(), SolutionMoistureError> {
        self.original.admit(cx)
    }
    fn advance(
        &self,
        dt_s: f64,
        _temperature: f64,
        activity: f64,
        cx: &Cx<'_>,
    ) -> Result<(Self::Field, Self::Report), SolutionMoistureError> {
        self.original.finish_stage(dt_s, self.known, activity, cx)
    }
    fn water_change(&self, next: &Self::Field) -> f64 {
        self.original.water_change(next)
    }
    fn water_mass_kg(&self, endpoint: Option<&Self::Field>) -> f64 {
        self.original.water_mass_kg(endpoint)
    }
}

impl FiniteMoistureSolution {
    /// Couple the original heterogeneous slab and finite solution with SDIRK2.
    /// Each stage evaluates actual bath activity and original nonlinear storage,
    /// sorption, diffusivity, and exterior/interface resistance together.
    ///
    /// The endpoint report includes both RK-weighted face transfers once; the
    /// first-stage report is diagnostic. Sources, temperature, reference
    /// geometry and solute totals remain original. Second-order accuracy assumes
    /// a smooth source/phase branch; kinetic crystal evolution is refused.
    /// The existing equal-pressure vapor bridge does not imply immersion,
    /// sorption heat or osmotic traction. Failure publishes neither inventory.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_layered_sdirk2(
        &self,
        moisture: &MoistureLayeredSlab,
        dt_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<LayeredSolutionSdirk2Candidate, SolutionMoistureError> {
        single(
            self.advance_shared_layered_sdirk2(
                &[LayeredMoistureParticipant {
                    id: 0,
                    moisture,
                    exposed,
                    transport,
                    nonlinear,
                }],
                dt_s,
                moisture.activity().len(),
                controls,
                cx,
            )?,
            cx,
        )
    }

    /// Share both finite-solution stages among complete original layered slabs.
    /// Canonical specimen IDs, local source/material inventories and one total
    /// solvent balance retain opposing uptake and release without extra debits.
    pub fn advance_shared_layered_sdirk2(
        &self,
        participants: &[LayeredMoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedLayeredSolutionSdirk2Candidate, SolutionMoistureError> {
        self.advance_layered_sdirk_participants(
            participants,
            dt_s,
            maximum_total_cells,
            controls,
            cx,
        )
    }

    /// Couple a coated cylinder/sphere and finite solution with native SDIRK2.
    /// Exact curved paths, unequal shell masses and each material's sorption
    /// remain in the original layered owner. The center of a solid body stays
    /// sealed, and interior material joins never become reservoir ports.
    ///
    /// The accuracy, stage accounting and fixed thermodynamic/geometry scope
    /// of [`Self::advance_layered_sdirk2`] apply unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_layered_radial_sdirk2(
        &self,
        moisture: &MoistureLayeredRadial,
        dt_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<LayeredRadialSolutionSdirk2Candidate, SolutionMoistureError> {
        single(
            self.advance_shared_layered_radial_sdirk2(
                &[LayeredRadialMoistureParticipant {
                    id: 0,
                    moisture,
                    exposed,
                    transport,
                    nonlinear,
                }],
                dt_s,
                moisture.activity().len(),
                controls,
                cx,
            )?,
            cx,
        )
    }

    /// Original coated radial bodies sharing one finite solution at both stages.
    /// The final root always uses original bath/storage, not the intermediate
    /// solution. All native source, phase, clock and aggregate cell limits apply.
    pub fn advance_shared_layered_radial_sdirk2(
        &self,
        participants: &[LayeredRadialMoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedLayeredRadialSolutionSdirk2Candidate, SolutionMoistureError> {
        self.advance_layered_sdirk_participants(
            participants,
            dt_s,
            maximum_total_cells,
            controls,
            cx,
        )
    }

    fn advance_layered_sdirk_participants<P: LayeredSdirkParticipant>(
        &self,
        participants: &[P],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: SolutionMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<
        SharedSolutionMoistureCandidate<P::Field, MoistureLayeredSdirk2Report>,
        SolutionMoistureError,
    > {
        poll(cx)?;
        crate::moisture::solution::shared::admit_count(participants.len(), maximum_total_cells)?;
        if self.kinetic_solubility().is_some() {
            return Err(SolutionMoistureError::Invalid(
                "second-order solution coupling requires algebraic phase partitioning",
            ));
        }
        let implicit_dt = kernel::sdirk_interval(participants[0].time_s(), dt_s)?;
        let (_, first_solids, first_report) = self
            .advance_participants(participants, implicit_dt, maximum_total_cells, controls, cx)?
            .into_parts();
        let mut ordered: Vec<_> = participants.iter().collect();
        ordered.sort_unstable_by_key(|body| body.id());
        let mut loads = Vec::with_capacity(ordered.len());
        for (body, first) in ordered.iter().zip(&first_solids) {
            poll(cx)?;
            debug_assert_eq!(body.id(), first.id);
            loads.push(body.stage_load(
                dt_s,
                &first.moisture,
                first_report.water_activity[1],
                cx,
            )?);
        }
        // Retain diagnostics, not complete intermediate material bodies,
        // while the final root allocates its own native candidates.
        let first_reports: Vec<_> = first_solids
            .into_iter()
            .map(|solid| (solid.id, solid.transport))
            .collect();
        let staged: Vec<_> = ordered
            .iter()
            .zip(&loads)
            .map(|(original, known)| FinalLayeredStage {
                original: *original,
                known,
            })
            .collect();
        // Original bath and body storage are deliberate. Only retained native
        // face loads carry stage-one evolution into this endpoint root.
        let (reservoir, solids, mut step) = self
            .advance_participants(&staged, dt_s, maximum_total_cells, controls, cx)?
            .into_parts();
        step.iterations = step.iterations.checked_add(first_report.iterations).ok_or(
            SolutionMoistureError::Invalid("layered SDIRK root work count overflow"),
        )?;
        let solids = solids
            .into_iter()
            .zip(first_reports)
            .map(|(last, (first_id, first_report))| {
                debug_assert_eq!(last.id, first_id);
                SharedMoistureSolid {
                    id: last.id,
                    moisture: last.moisture,
                    transport: MoistureLayeredSdirk2Report {
                        endpoint: last.transport,
                        first_stage: first_report,
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

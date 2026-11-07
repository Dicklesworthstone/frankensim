//! Both native layered diffusion stages coupled to one original finite gas.
//!
//! First-stage humidity supplies the physical face loads. The final root uses
//! ORIGINAL gas and material storage, so every weighted water transfer and the
//! whole-interval gas energy are counted once. All source, carrier and total
//! mass gates remain in the original atmosphere owner.

use super::*;
use crate::moisture::{
    layered::{
        MoistureLayeredSdirk2Report,
        kernel::{self, StageLoad},
    },
    radial::layered::{
        MoistureLayeredRadial,
        atmosphere::{self as radial_air, LayeredRadialMoistureParticipant},
    },
};

pub type LayeredAtmosphereSdirk2Step = AtmosphereMoistureStep<MoistureLayeredSdirk2Report>;
pub type LayeredAtmosphereSdirk2Candidate =
    AtmosphereMoistureCandidate<MoistureLayeredSlab, MoistureLayeredSdirk2Report>;
pub type SharedLayeredAtmosphereSdirk2Candidate =
    SharedAtmosphereMoistureCandidate<MoistureLayeredSlab, MoistureLayeredSdirk2Report>;
pub type LayeredRadialAtmosphereSdirk2Candidate =
    AtmosphereMoistureCandidate<MoistureLayeredRadial, MoistureLayeredSdirk2Report>;
pub type SharedLayeredRadialAtmosphereSdirk2Candidate =
    SharedAtmosphereMoistureCandidate<MoistureLayeredRadial, MoistureLayeredSdirk2Report>;

/// Only native geometry calls differ. The existing shared BE entry point owns
/// all ID, clock, aggregate-cell, source/environment and reservoir admission.
trait LayeredGasStage: Sized {
    type Field;
    fn id(&self) -> u64;
    fn time_s(&self) -> f64;
    fn sealed(&self) -> bool;
    fn first(
        atmosphere: &FiniteMoistureAtmosphere,
        participants: &[Self],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<
        SharedAtmosphereMoistureCandidate<Self::Field, MoistureLayeredReport>,
        AtmosphereMoistureError,
    >;
    fn load(
        &self,
        dt_s: f64,
        first: &Self::Field,
        rh: f64,
        cx: &Cx<'_>,
    ) -> Result<StageLoad, MoistureError>;
    fn finish(
        &self,
        dt_s: f64,
        known: &StageLoad,
        rh: f64,
        cx: &Cx<'_>,
    ) -> Result<(Self::Field, MoistureLayeredReport), MoistureError>;
    fn water_change(&self, next: &Self::Field) -> f64;
}

fn exposed(exposure: [bool; 2], transfer: &[MoistureActivityTransfer; 2]) -> [bool; 2] {
    std::array::from_fn(|side| exposure[side] && transfer[side] != MoistureActivityTransfer::Sealed)
}
fn boundary(exposure: [bool; 2], rh: f64) -> [MoistureLayeredFace; 2] {
    exposure.map(|face| {
        if face {
            MoistureLayeredFace::PrescribedRelativeHumidity(rh)
        } else {
            MoistureLayeredFace::Sealed
        }
    })
}

impl LayeredGasStage for LayeredMoistureParticipant<'_> {
    type Field = MoistureLayeredSlab;
    fn id(&self) -> u64 {
        self.id
    }
    fn time_s(&self) -> f64 {
        self.moisture.time_s()
    }
    fn sealed(&self) -> bool {
        exposed(self.exposed, self.moisture.exterior_transfer()) == [false; 2]
    }
    fn first(
        atmosphere: &FiniteMoistureAtmosphere,
        participants: &[Self],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedLayeredAtmosphereMoistureCandidate, AtmosphereMoistureError> {
        atmosphere.advance_shared_layered(participants, dt_s, maximum_total_cells, controls, cx)
    }
    fn load(
        &self,
        dt_s: f64,
        first: &Self::Field,
        rh: f64,
        cx: &Cx<'_>,
    ) -> Result<StageLoad, MoistureError> {
        let faces = boundary(exposed(self.exposed, self.moisture.exterior_transfer()), rh);
        self.moisture.sdirk_load(dt_s, first, &faces, cx)
    }
    fn finish(
        &self,
        dt_s: f64,
        known: &StageLoad,
        rh: f64,
        cx: &Cx<'_>,
    ) -> Result<(Self::Field, MoistureLayeredReport), MoistureError> {
        let faces = boundary(exposed(self.exposed, self.moisture.exterior_transfer()), rh);
        self.moisture
            .advance_sdirk_stage(dt_s, &faces, known, self.transport, self.nonlinear, cx)
    }
    fn water_change(&self, next: &Self::Field) -> f64 {
        super::water_change(self.moisture, next)
    }
}

impl LayeredGasStage for LayeredRadialMoistureParticipant<'_> {
    type Field = MoistureLayeredRadial;
    fn id(&self) -> u64 {
        self.id
    }
    fn time_s(&self) -> f64 {
        self.moisture.time_s()
    }
    fn sealed(&self) -> bool {
        exposed(self.exposed, self.moisture.exterior_transfer()) == [false; 2]
    }
    fn first(
        atmosphere: &FiniteMoistureAtmosphere,
        participants: &[Self],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<radial_air::SharedLayeredRadialAtmosphereMoistureCandidate, AtmosphereMoistureError>
    {
        atmosphere.advance_shared_layered_radial(
            participants,
            dt_s,
            maximum_total_cells,
            controls,
            cx,
        )
    }
    fn load(
        &self,
        dt_s: f64,
        first: &Self::Field,
        rh: f64,
        cx: &Cx<'_>,
    ) -> Result<StageLoad, MoistureError> {
        let faces = boundary(exposed(self.exposed, self.moisture.exterior_transfer()), rh);
        self.moisture.sdirk_load(dt_s, first, &faces, cx)
    }
    fn finish(
        &self,
        dt_s: f64,
        known: &StageLoad,
        rh: f64,
        cx: &Cx<'_>,
    ) -> Result<(Self::Field, MoistureLayeredReport), MoistureError> {
        let faces = boundary(exposed(self.exposed, self.moisture.exterior_transfer()), rh);
        self.moisture
            .advance_sdirk_stage(dt_s, &faces, known, self.transport, self.nonlinear, cx)
    }
    fn water_change(&self, next: &Self::Field) -> f64 {
        radial_air::water_change(self.moisture, next)
    }
}

fn single<F>(
    candidate: SharedAtmosphereMoistureCandidate<F, MoistureLayeredSdirk2Report>,
    surface: impl FnOnce(
        &MoistureLayeredReport,
    ) -> Result<[Option<ResolvedMaterialStatePoint>; 2], MoistureError>,
    cx: &Cx<'_>,
) -> Result<AtmosphereMoistureCandidate<F, MoistureLayeredSdirk2Report>, AtmosphereMoistureError> {
    let (atmosphere, mut solids, shared) = candidate.into_parts();
    let solid = solids
        .pop()
        .expect("one admitted layered atmosphere participant");
    let sorption = surface(&solid.transport.endpoint)?;
    let step = AtmosphereMoistureStep {
        sorption,
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

impl FiniteMoistureAtmosphere {
    /// Couple second-order native heterogeneous slab diffusion to finite vapor.
    ///
    /// Both stages resolve their actual gas humidity and each original material's
    /// nonlinear sorption, diffusivity and interface/exterior resistance. The
    /// endpoint receipt contains the full RK-weighted transfers once; the first
    /// transport stage is diagnostic and must not be debited again.
    ///
    /// Temperature, gas volume and body geometry stay fixed. The original gas
    /// source, carrier/total-mass and whole-interval caloric checks remain active.
    /// This does not add condensation, sorption heat or a temporal certificate.
    /// Inadmissible stages refuse atomically; callers may shorten the interval.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_layered_sdirk2(
        &self,
        moisture: &MoistureLayeredSlab,
        dt_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<LayeredAtmosphereSdirk2Candidate, AtmosphereMoistureError> {
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
            |report| super::surface_sorption(moisture, report, cx),
            cx,
        )
    }

    /// One common humidity per stage for independently sourced layered slabs.
    /// Canonical IDs, original capacities and a single vapor balance preserve
    /// opposing uptake/release. Every stage starts from original accepted storage.
    pub fn advance_shared_layered_sdirk2(
        &self,
        participants: &[LayeredMoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedLayeredAtmosphereSdirk2Candidate, AtmosphereMoistureError> {
        self.advance_layered_gas_stages(participants, dt_s, maximum_total_cells, controls, cx)
    }

    /// Native two-stage transport for a coated cylinder/sphere and finite gas.
    /// Each layer retains its source and actual shell masses; the stage load
    /// retains exact curved paths and physical material-interface transfers.
    /// A solid center stays sealed and interior joins create no gas ports.
    ///
    /// The complete-stage accounting and fixed physical scope of
    /// [Self::advance_layered_sdirk2] apply unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_layered_radial_sdirk2(
        &self,
        moisture: &MoistureLayeredRadial,
        dt_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<LayeredRadialAtmosphereSdirk2Candidate, AtmosphereMoistureError> {
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
            |report| radial_air::surface_sorption(moisture, report, cx),
            cx,
        )
    }

    /// Independently sourced coated bodies sharing both native finite-gas stages.
    /// Actual unequal shell storage weights the compensated water account.
    /// The final gas-energy increment is original-to-final, not a sum that
    /// double-counts the temporary first-stage atmosphere.
    pub fn advance_shared_layered_radial_sdirk2(
        &self,
        participants: &[LayeredRadialMoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedLayeredRadialAtmosphereSdirk2Candidate, AtmosphereMoistureError> {
        self.advance_layered_gas_stages(participants, dt_s, maximum_total_cells, controls, cx)
    }

    fn advance_layered_gas_stages<P: LayeredGasStage>(
        &self,
        participants: &[P],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<
        SharedAtmosphereMoistureCandidate<P::Field, MoistureLayeredSdirk2Report>,
        AtmosphereMoistureError,
    > {
        poll(cx)?;
        admit_participants(participants.len(), maximum_total_cells)?;
        let start = participants[0].time_s();
        let implicit_dt = kernel::sdirk_interval(start, dt_s)?;
        let (_, first_solids, first_report) = P::first(
            self,
            participants,
            implicit_dt,
            maximum_total_cells,
            controls,
            cx,
        )?
        .into_parts();
        let mut ordered: Vec<_> = participants.iter().collect();
        ordered.sort_unstable_by_key(|body| body.id());
        let mut loads = Vec::with_capacity(ordered.len());
        for (body, first) in ordered.iter().zip(&first_solids) {
            poll(cx)?;
            debug_assert_eq!(body.id(), first.id);
            loads.push(body.load(dt_s, &first.moisture, first_report.relative_humidity[1], cx)?);
        }
        // Drop intermediate bodies before the final root allocates candidates.
        // Only native weighted face loads and first-stage diagnostics survive.
        let first_reports: Vec<_> = first_solids
            .into_iter()
            .map(|solid| (solid.id, solid.transport))
            .collect();
        let all_sealed = ordered.iter().all(|body| body.sealed());
        let (atmosphere, last_solids, mut step) =
            self.advance_exchange(start, all_sealed, controls, cx, |rh| {
                let mut solids = Vec::with_capacity(ordered.len());
                for (body, known) in ordered.iter().zip(&loads) {
                    poll(cx)?;
                    let (moisture, transport) = body.finish(dt_s, known, rh, cx)?;
                    solids.push(SharedMoistureSolid {
                        id: body.id(),
                        moisture,
                        transport,
                    });
                }
                let change = water_sum(
                    ordered
                        .iter()
                        .zip(&solids)
                        .map(|(old, next)| old.water_change(&next.moisture)),
                );
                let end_time = solids[0].transport.interval_s[1];
                Ok((solids, change, end_time))
            })?;
        step.iterations = step.iterations.checked_add(first_report.iterations).ok_or(
            AtmosphereMoistureError::Invalid("layered SDIRK root work count overflow"),
        )?;
        step.maximum_carrier_mass_residual_kg = step
            .maximum_carrier_mass_residual_kg
            .max(first_report.maximum_carrier_mass_residual_kg);
        let solids = last_solids
            .into_iter()
            .zip(first_reports)
            .map(|(last, (first_id, first_stage))| {
                debug_assert_eq!(last.id, first_id);
                SharedMoistureSolid {
                    id: last.id,
                    moisture: last.moisture,
                    transport: MoistureLayeredSdirk2Report {
                        endpoint: last.transport,
                        first_stage,
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

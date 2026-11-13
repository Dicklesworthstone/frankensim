//! Radial diffusion against original finite, well-mixed vapor inventories.
//!
//! Every humidity trial restarts the accepted radial field. Actual shell dry
//! masses weight its water increments; neither a slab nor a uniform-cell
//! inventory substitutes for the curved geometry. The existing atmosphere
//! owner reconstructs gas composition and pressure and checks water, carrier
//! and total mass before publishing the gas and solid endpoints together.
//!
//! Shared chambers debit their inventory once for all specimens. Paired
//! chambers retain the signed inner/outer face accounts, including through-flow
//! with no net solid uptake. Temperature, gas volumes and reference geometry
//! remain fixed. Gas internal-energy change is not total sorption heat; carrier
//! permeation, condensation and thermal feedback remain outside this model.

use super::MoistureRadial;
use crate::moisture::{
    MoistureNonlinearConfig, MoistureStepConfig, MoistureStepReport,
    atmosphere::{
        AtmosphereMoistureCandidate, AtmosphereMoistureError, AtmosphereMoistureSolve,
        AtmosphereMoistureStep, DryBasisSorption, FiniteMoistureAtmosphere,
        SharedAtmosphereMoistureCandidate, SharedMoistureSolid, admit_participants,
        pair::{
            AtmosphereMoistureFace, MoistureAtmospherePort, PairedAtmosphereMoistureCandidate,
            PairedAtmosphereMoistureSolve, advance_pair_with,
        },
    },
    poll,
    reservoir::water_sum,
    surface::ResolvedMoistureFace,
};
use fs_exec::Cx;

mod sdirk;
pub use sdirk::{
    RadialAtmosphereSdirk2Candidate, RadialAtmosphereSdirk2Step,
    SharedRadialAtmosphereSdirk2Candidate,
};

/// One original gas inventory and radial field accepted at a common endpoint.
pub type RadialAtmosphereMoistureCandidate =
    AtmosphereMoistureCandidate<MoistureRadial, MoistureStepReport>;
/// One common finite atmosphere and its independently sourced radial specimens.
pub type SharedRadialAtmosphereMoistureCandidate =
    SharedAtmosphereMoistureCandidate<MoistureRadial, MoistureStepReport>;
/// Two independent chambers and the original radial field accepted together.
pub type PairedRadialAtmosphereMoistureCandidate =
    PairedAtmosphereMoistureCandidate<MoistureRadial, MoistureStepReport>;

/// A physical cylinder or sphere sharing one finite atmospheric inventory.
pub struct RadialMoistureParticipant<'a> {
    /// Stable, distinct physical-inventory ID; determines summation order.
    pub id: u64,
    /// Original accepted radial transport field.
    pub moisture: &'a MoistureRadial,
    /// Inner/outer surface laws; None is sealed. A solid center must be None.
    pub faces: &'a [Option<DryBasisSorption>; 2],
    /// Native diffusion work, residual and content-increment limits.
    pub transport: MoistureStepConfig,
    /// Native endpoint-dependent Kirchhoff iteration controls.
    pub nonlinear: MoistureNonlinearConfig,
}

fn water_change(old: &MoistureRadial, next: &MoistureRadial) -> f64 {
    // Subtract represented contents before scaling: subtracting large absolute
    // totals can erase small exchanges between unequal radial shell inventories.
    water_sum(
        old.moisture()
            .iter()
            .zip(next.moisture())
            .zip(old.dry_masses_kg())
            .map(|((&before, &after), &dry)| (after - before) * dry),
    )
}

fn admit(
    moisture: &MoistureRadial,
    temperature: f64,
    exposed: [bool; 2],
) -> Result<(), AtmosphereMoistureError> {
    if moisture.temperature_kelvin().to_bits() != temperature.to_bits() {
        return Err(AtmosphereMoistureError::Invalid(
            "radial specimen and finite atmospheres need the same prescribed temperature",
        ));
    }
    if moisture.geometry().faces_m[0] == 0.0 && exposed[0] {
        return Err(AtmosphereMoistureError::Invalid(
            "a solid radial center has zero area and must remain sealed",
        ));
    }
    Ok(())
}

impl FiniteMoistureAtmosphere {
    /// Couple a cylinder or sphere to this one finite vapor inventory. The
    /// actual endpoint RH drives the original sorption and finite surface laws;
    /// interior transport retains its nonlinear radial Kirchhoff resistance.
    /// Every rejected source, transport or gas trial leaves both inputs intact.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_radial(
        &self,
        moisture: &MoistureRadial,
        dt_s: f64,
        faces: &[Option<DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<RadialAtmosphereMoistureCandidate, AtmosphereMoistureError> {
        let participant = RadialMoistureParticipant {
            id: 0,
            moisture,
            faces,
            transport,
            nonlinear,
        };
        let (atmosphere, mut solids, shared) = self
            .advance_shared_radial(
                &[participant],
                dt_s,
                moisture.moisture().len(),
                controls,
                cx,
            )?
            .into_parts();
        let solid = solids.pop().expect("one admitted radial specimen");
        let step = AtmosphereMoistureStep {
            sorption: solid.transport.surface_sorption.clone(),
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

    /// Simultaneous uptake and release by any admitted radial specimens share
    /// one endpoint RH and one gas debit. Stable IDs fix ordering; accepted
    /// clocks and temperatures must agree. The total cell budget includes sealed
    /// specimens, whose internal diffusion still advances during the interval.
    pub fn advance_shared_radial(
        &self,
        participants: &[RadialMoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedRadialAtmosphereMoistureCandidate, AtmosphereMoistureError> {
        poll(cx)?;
        admit_participants(participants.len(), maximum_total_cells)?;
        let mut ordered: Vec<_> = participants.iter().collect();
        ordered.sort_unstable_by_key(|body| body.id);
        if ordered.windows(2).any(|pair| pair[0].id == pair[1].id) {
            return Err(AtmosphereMoistureError::Invalid(
                "shared-atmosphere radial specimen IDs must be unique",
            ));
        }
        let time_s = ordered[0].moisture.time_s();
        let temperature = self.gas().state().state().temperature;
        let mut cells = 0usize;
        for body in &ordered {
            poll(cx)?;
            cells = cells
                .checked_add(body.moisture.moisture().len())
                .filter(|&n| n <= maximum_total_cells)
                .ok_or(AtmosphereMoistureError::Invalid(
                    "shared-atmosphere radial cell budget exhausted",
                ))?;
            if body.moisture.time_s().to_bits() != time_s.to_bits() {
                return Err(AtmosphereMoistureError::Invalid(
                    "radial specimens need the same accepted clock",
                ));
            }
            admit(
                body.moisture,
                temperature,
                body.faces.each_ref().map(Option::is_some),
            )?;
        }
        let all_sealed = ordered
            .iter()
            .all(|body| body.faces.iter().all(Option::is_none));
        let (atmosphere, solids, step) =
            self.advance_exchange(time_s, all_sealed, controls, cx, |rh| {
                let mut solids = Vec::with_capacity(ordered.len());
                for body in &ordered {
                    poll(cx)?;
                    let mut boundary = std::array::from_fn(|_| ResolvedMoistureFace::default());
                    for (side, law) in body.faces.iter().enumerate() {
                        if let Some(law) = law {
                            boundary[side] = law.resolve(temperature, rh)?;
                        }
                    }
                    let (moisture, transport) = body.moisture.advance_stage(
                        dt_s,
                        dt_s,
                        boundary,
                        body.transport,
                        body.nonlinear,
                        None,
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
                        .map(|(old, next)| water_change(old.moisture, &next.moisture)),
                );
                let end_time = solids[0].moisture.time_s();
                Ok((solids, change, end_time))
            })?;
        Ok(SharedAtmosphereMoistureCandidate::from_accepted(
            atmosphere, solids, step,
        ))
    }

    /// Exchange through an annulus or spherical shell between two independent
    /// finite chambers. The receiver is First; `second` is Second. Each radial
    /// face names its one chamber, and its signed native flux debits only that
    /// chamber even when opposing transfers give zero net solid uptake.
    ///
    /// The original nested RH solver bounds all field evaluations and may
    /// refuse unbracketed or unsupported probe paths. No interval, tolerance or
    /// source domain is widened, and all three accepted inputs remain unchanged
    /// until individual and whole-system water/mass checks succeed.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_paired_radial(
        &self,
        second: &Self,
        moisture: &MoistureRadial,
        dt_s: f64,
        faces: &[Option<AtmosphereMoistureFace>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: PairedAtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<PairedRadialAtmosphereMoistureCandidate, AtmosphereMoistureError> {
        poll(cx)?;
        let temperature = self.gas().state().state().temperature;
        admit(moisture, temperature, faces.each_ref().map(Option::is_some))?;
        let face_atmospheres = faces
            .each_ref()
            .map(|face| face.as_ref().map(|f| f.atmosphere));
        advance_pair_with(
            [self, second],
            moisture.time_s(),
            face_atmospheres,
            controls,
            cx,
            |rh| {
                let mut boundary = std::array::from_fn(|_| ResolvedMoistureFace::default());
                for (side, face) in faces.iter().enumerate() {
                    if let Some(face) = face {
                        boundary[side] = face.sorption.resolve(
                            temperature,
                            match face.atmosphere {
                                MoistureAtmospherePort::First => rh[0],
                                MoistureAtmospherePort::Second => rh[1],
                            },
                        )?;
                    }
                }
                let (next, report) =
                    moisture.advance_stage(dt_s, dt_s, boundary, transport, nonlinear, None, cx)?;
                let transfers = report.boundary_water_kg;
                let end_time = next.time_s();
                Ok((next, report, transfers, end_time))
            },
            |next| water_change(moisture, next),
        )
    }
}

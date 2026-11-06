//! Coupled two-stage diffusion against one evolving finite vapor inventory.
//!
//! Each stage solves all exposed fields and one common RH together. Stage two
//! keeps the ORIGINAL solid and gas storage references and adds the actual
//! first-stage flux load. The scalar water equation therefore includes the
//! explicit RK transfer exactly once, without a second reservoir debit.

use super::*;
use crate::moisture::cartesian::MoistureCartesianSdirk2Report;
use crate::moisture::sdirk::{GAMMA, StageLoad};
use std::borrow::Cow;

/// Original-to-final finite-vapor account with both spatial stage receipts.
pub type CartesianAtmosphereSdirk2Step = AtmosphereMoistureStep<MoistureCartesianSdirk2Report, 6>;
/// One accepted atmosphere and Cartesian field after both coupled stages.
pub type CartesianAtmosphereSdirk2Candidate =
    AtmosphereMoistureCandidate<MoistureCartesian, MoistureCartesianSdirk2Report, 6>;
/// Multiple original Cartesian fields sharing one second-order vapor account.
pub type SharedCartesianAtmosphereSdirk2Candidate =
    SharedAtmosphereMoistureCandidate<MoistureCartesian, MoistureCartesianSdirk2Report>;

// Borrow the original field/load during scalar trials. Only a successful native
// solve owns a proposed field; neither storage nor the load is copied per trial.
struct CartesianStage<'a> {
    field: Cow<'a, MoistureCartesian>,
    known: &'a StageLoad<6>,
    implicit_dt_s: f64,
}
impl ReservoirField<6> for CartesianStage<'_> {
    type Report = MoistureCartesianReport;
    type Controls = (MoistureStepConfig, MoistureNonlinearConfig);

    fn inventory(&self) -> super::super::binding::MoistureField<'_> {
        self.field.field()
    }

    fn advance_resolved(
        &self,
        dt_s: f64,
        boundaries: [ResolvedMoistureFace; 6],
        controls: Self::Controls,
        cx: &Cx<'_>,
    ) -> Result<(Self, Self::Report), MoistureError> {
        let (field, report) = self.field.advance_stage(
            dt_s,
            self.implicit_dt_s,
            boundaries,
            controls.0,
            controls.1,
            Some(self.known),
            cx,
        )?;
        Ok((
            Self {
                field: Cow::Owned(field),
                known: self.known,
                implicit_dt_s: self.implicit_dt_s,
            },
            report,
        ))
    }
}

impl FiniteMoistureAtmosphere {
    /// Second-order Cartesian diffusion coupled to the actual finite vapor.
    ///
    /// Sorption, surface emission, nonlinear principal diffusivity and gas
    /// composition are resolved at each stage's own RH. The endpoint transport
    /// report contains RK-weighted physical transfers; its first-stage report
    /// is diagnostic and must not be added again. Outer iterations count both
    /// scalar solves; the gas and water accounts cover only the full interval.
    ///
    /// Temperature and volume remain prescribed, with the original isothermal
    /// caloric scope. No condensation, sorption heat or global error bound is
    /// introduced. Inadmissible stages or source-bracket trials refuse the
    /// entire proposal and can require a shorter interval.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_cartesian_sdirk2(
        &self,
        moisture: &MoistureCartesian,
        dt_s: f64,
        faces: &[Option<DryBasisSorption>; 6],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<CartesianAtmosphereSdirk2Candidate, AtmosphereMoistureError> {
        let participant = CartesianMoistureParticipant {
            id: 0,
            moisture,
            faces,
            transport,
            nonlinear,
        };
        let (atmosphere, mut solids, shared) = self
            .advance_shared_cartesian_sdirk2(
                &[participant],
                dt_s,
                moisture.moisture().len(),
                controls,
                cx,
            )?
            .into_parts();
        let solid = solids.pop().expect("one admitted Cartesian participant");
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
        Ok(AtmosphereMoistureCandidate {
            atmosphere,
            moisture: solid.moisture,
            step,
        })
    }

    /// Advance heterogeneous 3D specimens with one shared RH at each SDIRK stage.
    ///
    /// Each specimen retains its own six source-backed faces and anisotropic
    /// nonlinear diffusion owner. Stable-ID ordering, common accepted clocks,
    /// temperature, aggregate cell budgets and native gas/species checks are
    /// the same as ordinary shared finite-vapor evolution. Simultaneous uptake
    /// and desorption share one compensated water balance.
    ///
    /// Stage-one storage is diagnostic. Both stage equations reference the
    /// original accepted inventories, so no specimen can spend the vapor twice.
    pub fn advance_shared_cartesian_sdirk2(
        &self,
        participants: &[CartesianMoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedCartesianAtmosphereSdirk2Candidate, AtmosphereMoistureError> {
        poll(cx)?;
        admit_participants(participants.len(), maximum_total_cells)?;
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
                "finite-vapor SDIRK needs two representable forward stage times",
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
        let temperature = self.gas.state().state().temperature;
        let first_rh = first_report.relative_humidity[1];
        let mut loads = Vec::with_capacity(ordered.len());
        for (body, first) in ordered.iter().zip(&first_solids) {
            poll(cx)?;
            debug_assert_eq!(body.id, first.id);
            let mut boundaries = std::array::from_fn(|_| ResolvedMoistureFace::default());
            for (side, law) in body.faces.iter().enumerate() {
                if let Some(law) = law {
                    boundaries[side] = law.resolve(temperature, first_rh)?;
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
            .map(|(body, known)| CartesianStage {
                field: Cow::Borrowed(body.moisture),
                known,
                implicit_dt_s,
            })
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

        // Intentionally solve from self, not the diagnostic stage-one gas.
        let (atmosphere, solids, mut step) = self
            .advance_fields(&staged, dt_s, maximum_total_cells, controls, cx)?
            .into_parts();
        step.iterations += first_report.iterations;
        step.maximum_carrier_mass_residual_kg = step
            .maximum_carrier_mass_residual_kg
            .max(first_report.maximum_carrier_mass_residual_kg);
        let solids = solids
            .into_iter()
            .zip(first_solids)
            .map(|(last, first)| {
                debug_assert_eq!(last.id, first.id);
                SharedMoistureSolid {
                    id: last.id,
                    moisture: last.moisture.field.into_owned(),
                    transport: MoistureCartesianSdirk2Report {
                        endpoint: last.transport,
                        first_stage: first.transport,
                    },
                }
            })
            .collect();
        poll(cx)?;
        Ok(SharedAtmosphereMoistureCandidate {
            atmosphere,
            solids,
            step,
        })
    }
}

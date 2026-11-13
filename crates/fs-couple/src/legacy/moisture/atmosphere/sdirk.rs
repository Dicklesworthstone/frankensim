//! Coupled two-stage diffusion against one evolving finite vapor inventory.
//!
//! Each stage solves all exposed fields and one common RH together. Stage two
//! keeps the ORIGINAL solid and gas storage references and adds the actual
//! first-stage flux load. The scalar water equation therefore includes the
//! explicit RK transfer exactly once, without a second reservoir debit.

use super::*;
use crate::moisture::cartesian::MoistureCartesianSdirk2Report;
use crate::moisture::reservoir::{CartesianSdirkStage, SlabSdirkStage};
use crate::moisture::sdirk::{GAMMA, MoistureSdirk2Report};

/// Original-to-final finite-vapor account with both native slab stage reports.
pub type SlabAtmosphereSdirk2Step = AtmosphereMoistureStep<MoistureSdirk2Report>;
/// One accepted finite atmosphere and conservative tridiagonal slab endpoint.
pub type SlabAtmosphereSdirk2Candidate =
    AtmosphereMoistureCandidate<MoistureSlab, MoistureSdirk2Report>;
/// Multiple native slabs accepted with one shared second-order vapor account.
pub type SharedSlabAtmosphereSdirk2Candidate =
    SharedAtmosphereMoistureCandidate<MoistureSlab, MoistureSdirk2Report>;

/// Original-to-final finite-vapor account with both spatial stage receipts.
pub type CartesianAtmosphereSdirk2Step = AtmosphereMoistureStep<MoistureCartesianSdirk2Report, 6>;
/// One accepted atmosphere and Cartesian field after both coupled stages.
pub type CartesianAtmosphereSdirk2Candidate =
    AtmosphereMoistureCandidate<MoistureCartesian, MoistureCartesianSdirk2Report, 6>;
/// Multiple original Cartesian fields sharing one second-order vapor account.
pub type SharedCartesianAtmosphereSdirk2Candidate =
    SharedAtmosphereMoistureCandidate<MoistureCartesian, MoistureCartesianSdirk2Report>;

impl FiniteMoistureAtmosphere {
    /// Second-order moisture diffusion on the native one-dimensional slab,
    /// coupled to this atmosphere's actual finite water inventory at both stages.
    ///
    /// Each source-dependent Picard iterate retains the O(N) tridiagonal solve.
    /// Sorption, surface emission and D(M) are evaluated at each stage's own
    /// accepted RH and cell moisture, including Kirchhoff transport when selected.
    /// Both stage balances use the original water inventories; first-stage gas
    /// storage is diagnostic and is never debited twice. The returned slab keeps
    /// its original choice of ordinary time integrator.
    ///
    /// Temperature, gas volume and reference skeleton stay fixed. This does not
    /// add sorption heat, condensation or a global temporal-error certificate.
    /// Inadmissible stages/source trials refuse the whole proposal without
    /// clipping or first-order fallback; a shorter interval may be required.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_sdirk2(
        &self,
        slab: &MoistureSlab,
        dt_s: f64,
        faces: &[Option<DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SlabAtmosphereSdirk2Candidate, AtmosphereMoistureError> {
        let participant = MoistureParticipant {
            id: 0,
            slab,
            faces,
            transport,
        };
        let (atmosphere, mut solids, shared) = self
            .advance_shared_sdirk2(
                &[participant],
                dt_s,
                slab.moisture().len(),
                nonlinear,
                controls,
                cx,
            )?
            .into_parts();
        let solid = solids.pop().expect("one admitted slab participant");
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

    /// Advance native slabs through two stages with one common finite vapor
    /// inventory. Each specimen retains its own material, geometry, surface
    /// laws and original ordinary integration mode. The supplied nonlinear
    /// controls apply to both implicit stages of every specimen.
    ///
    /// Stable-ID ordering and compensated signed transfers preserve simultaneous
    /// uptake and drying, including exchanges with little net vapor change.
    /// Stage two uses ORIGINAL gas/solid storage and the weighted first-stage
    /// face loads once. The intermediate gas and stage reports are diagnostic.
    /// Source, clock, temperature, aggregate-cell and cancellation checks admit
    /// every endpoint together; an invalid specimen publishes no partial state.
    ///
    /// As for one slab, temperature and gas volume stay prescribed. SDIRK
    /// stability is not unconditional monotonicity; a shorter step may be
    /// required for source support or a physical second-stage moisture range.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_shared_sdirk2(
        &self,
        participants: &[MoistureParticipant<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        nonlinear: MoistureNonlinearConfig,
        controls: AtmosphereMoistureSolve,
        cx: &Cx<'_>,
    ) -> Result<SharedSlabAtmosphereSdirk2Candidate, AtmosphereMoistureError> {
        poll(cx)?;
        nonlinear.validate()?;
        admit_participants(participants.len(), maximum_total_cells)?;
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
                "finite-vapor slab SDIRK needs two representable forward stage times",
            )
            .into());
        }

        // Configure each original field once, not once per scalar RH trial.
        // Admit the aggregate allocation before copying any specimen state.
        let mut cells = 0usize;
        for body in &ordered {
            poll(cx)?;
            cells = cells
                .checked_add(body.slab.moisture().len())
                .filter(|count| *count <= maximum_total_cells)
                .ok_or_else(|| invalid("shared-atmosphere total cell budget exhausted"))?;
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

        let temperature = self.gas.state().state().temperature;
        let first_rh = first_report.relative_humidity[1];
        let mut loads = Vec::with_capacity(ordered.len());
        for ((body, field), first) in ordered.iter().zip(&staged).zip(&first_solids) {
            poll(cx)?;
            debug_assert_eq!(body.id, first.id);
            let mut boundaries = std::array::from_fn(|_| ResolvedMoistureFace::default());
            for (side, law) in body.faces.iter().enumerate() {
                if let Some(law) = law {
                    boundaries[side] = law.resolve(temperature, first_rh)?;
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
        // The accepted first-stage receipts and rates suffice from here.
        // Release its full diagnostic fields before the second scalar search.
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
        // Do not start from the diagnostic first-stage gas: that would spend
        // the same finite vapor a second time.
        let (atmosphere, solids, mut step) = self
            .advance_fields(&last_fields, dt_s, maximum_total_cells, controls, cx)?
            .into_parts();
        step.iterations += first_report.iterations;
        step.maximum_carrier_mass_residual_kg = step
            .maximum_carrier_mass_residual_kg
            .max(first_report.maximum_carrier_mass_residual_kg);
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
        Ok(SharedAtmosphereMoistureCandidate {
            atmosphere,
            solids,
            step,
        })
    }

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
                    moisture: last.moisture.into_field(),
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

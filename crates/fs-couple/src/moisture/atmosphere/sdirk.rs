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
        poll(cx)?;
        nonlinear.validate()?;
        let implicit_dt_s = GAMMA * dt_s;
        let first_time = slab.time_s() + implicit_dt_s;
        let end = slab.time_s() + dt_s;
        if !positive(dt_s) || !positive(implicit_dt_s) || !end.is_finite()
            || first_time <= slab.time_s() || first_time >= end
        {
            return Err(MoistureError::Invalid(
                "finite-vapor slab SDIRK needs two representable forward stage times",
            ).into());
        }
        let mut staged = slab.clone();
        staged.implicit_diffusivity = Some(nonlinear);
        if staged.kirchhoff_diffusivity.is_some() {
            staged.kirchhoff_diffusivity = Some(nonlinear);
        }
        let (_, first, first_report) = self.advance(
            &staged, implicit_dt_s, faces, transport, controls, cx,
        )?.into_parts();
        let mut boundaries = std::array::from_fn(|_| ResolvedMoistureFace::default());
        for (side, law) in faces.iter().enumerate() {
            poll(cx)?;
            if let Some(law) = law {
                boundaries[side] = law.resolve(
                    self.gas.state().state().temperature, first_report.relative_humidity[1],
                )?;
            }
        }
        let known = staged.sdirk_load(
            dt_s, &first, &first_report.transport, &boundaries, cx,
        )?;
        let field = SlabSdirkStage::new(&staged, &known, implicit_dt_s);
        let participant = ReservoirParticipant {
            id: 0, field: &field, faces, transport,
        };
        // Stage two retains ORIGINAL solid/gas storage plus the first-stage
        // face load; starting from the diagnostic first gas would double-debit.
        let (atmosphere, mut solids, shared) = self.advance_fields(
            &[participant], dt_s, slab.moisture().len(), controls, cx,
        )?.into_parts();
        let solid = solids.pop().expect("one admitted slab participant");
        let mut moisture = solid.moisture.into_field();
        moisture.implicit_diffusivity = slab.implicit_diffusivity;
        moisture.kirchhoff_diffusivity = slab.kirchhoff_diffusivity;
        let step = AtmosphereMoistureStep {
            sorption: solid.transport.surface_sorption.clone(),
            transport: MoistureSdirk2Report {
                endpoint: solid.transport,
                first_stage: first_report.transport,
            },
            relative_humidity: shared.relative_humidity,
            vapor_mass_kg: shared.vapor_mass_kg,
            water_to_solid_kg: shared.water_to_solids_kg,
            water_balance_residual_kg: shared.water_balance_residual_kg,
            maximum_carrier_mass_residual_kg: shared.maximum_carrier_mass_residual_kg
                .max(first_report.maximum_carrier_mass_residual_kg),
            gas_internal_energy_change_j: shared.gas_internal_energy_change_j,
            iterations: shared.iterations + first_report.iterations,
        };
        poll(cx)?;
        Ok(AtmosphereMoistureCandidate { atmosphere, moisture, step })
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

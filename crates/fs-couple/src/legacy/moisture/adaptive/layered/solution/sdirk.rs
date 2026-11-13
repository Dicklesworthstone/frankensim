//! Complete heterogeneous finite-bath durations with native coupled SDIRK2.

use super::*;

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
impl FiniteMoistureSolution {
    /// Complete layered slab diffusion and finite-solution feedback using the
    /// original second-order full-versus-two-half controller. Every native
    /// stage solves actual bath activity and each material's nonlinear storage,
    /// source diffusion and interface/surface resistance together.
    ///
    /// Local content/activity, actual solvent mass/activity and the L1 physical
    /// water discrepancy across all original layers independently enter the
    /// error comparison. Both implicit stages spend the work cap, including
    /// rejected probes. Only accepted fine endpoint transfers enter each layer
    /// and the exterior whole-body accounts; first-stage diagnostics are not
    /// debited again. One fixed original-to-final water tolerance applies.
    ///
    /// Geometry, temperature and exposed faces stay held. Second-order local
    /// estimates assume smooth source/phase branches and do not certify global
    /// accuracy. Kinetic phases retain the native SDIRK refusal. Any late
    /// failure or cancellation discards the complete bath/body duration.
    pub fn advance_layered_sdirk2_adaptive_to(
        &self,
        moisture: &MoistureLayeredSlab,
        target_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: SolutionMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            MoistureLayeredSlab,
            LayeredMoistureAdvanceReport<SolutionMoistureStep<MoistureLayeredSdirk2Report>>,
        ),
        MoistureAdvanceError<SolutionMoistureError>,
    > {
        checkpoint(cx)?;
        let initial = (self.clone(), moisture.clone());
        let ((reservoir, field), report) = evolve_fields_with_method(
            &initial,
            target_s,
            time,
            TimeMethod::Sdirk2,
            |state, dt| {
                let (bath, field, step) = state
                    .0
                    .advance_layered_sdirk2(
                        &state.1, dt, exposed, transport, nonlinear, physical, cx,
                    )?
                    .into_parts();
                Ok(((bath, field), LayeredStep(step)))
            },
            |state| {
                let mut snapshot = observe(&state.1);
                snapshot.reservoir_water_kg = Some(state.0.solvent_mass_kg());
                snapshot.activity = Some(state.0.solution().water_activity()?);
                Ok(snapshot)
            },
            |old, coarse, fine| {
                water_distance(
                    old.1
                        .layers()
                        .iter()
                        .zip(coarse.1.layers())
                        .zip(fine.1.layers())
                        .map(|((a, b), c)| (a.slab().field(), b.slab().field(), c.slab().field())),
                    time,
                    cx,
                )
            },
            cx,
        )?;
        Ok((
            reservoir,
            field,
            finish(report, time.water_balance_tolerance_kg, cx)?,
        ))
    }

    /// Complete finite-solution exchange with a coated sphere or cylinder using
    /// native second-order material and reservoir stages. Each shell retains
    /// its actual dry mass and curved resistance; every material keeps its own
    /// sorption, source support and interface partition.
    ///
    /// Error scales, full-duration conservation, stage work and failure
    /// semantics are those of [`Self::advance_layered_sdirk2_adaptive_to`]. A
    /// zero-area center cannot be exposed, and no kinetic phase step is relabeled
    /// as second order. The original vapor-equilibrium activity bridge remains.
    pub fn advance_layered_radial_sdirk2_adaptive_to(
        &self,
        moisture: &MoistureLayeredRadial,
        target_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: SolutionMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            MoistureLayeredRadial,
            LayeredMoistureAdvanceReport<SolutionMoistureStep<MoistureLayeredSdirk2Report>>,
        ),
        MoistureAdvanceError<SolutionMoistureError>,
    > {
        checkpoint(cx)?;
        let initial = (self.clone(), moisture.clone());
        let ((reservoir, field), report) = evolve_fields_with_method(
            &initial,
            target_s,
            time,
            TimeMethod::Sdirk2,
            |state, dt| {
                let (bath, field, step) = state
                    .0
                    .advance_layered_radial_sdirk2(
                        &state.1, dt, exposed, transport, nonlinear, physical, cx,
                    )?
                    .into_parts();
                Ok(((bath, field), LayeredStep(step)))
            },
            |state| {
                let mut snapshot = observe_radial(&state.1);
                snapshot.reservoir_water_kg = Some(state.0.solvent_mass_kg());
                snapshot.activity = Some(state.0.solution().water_activity()?);
                Ok(snapshot)
            },
            |old, coarse, fine| radial_water_distance(&old.1, &coarse.1, &fine.1, time, cx),
            cx,
        )?;
        Ok((
            reservoir,
            field,
            finish(report, time.water_balance_tolerance_kg, cx)?,
        ))
    }
}

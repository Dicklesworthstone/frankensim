//! The existing bounded duration controller on original layered SDIRK endpoints.

use super::*;

impl MoistureLayeredSlab {
    /// Evolve heterogeneous layers with their native second-order stages.
    /// Each local moisture and source activity enters the error comparison;
    /// `reservoir_water_absolute_tolerance_kg` additionally bounds the L1 sum
    /// of actual dry-mass-weighted cell discrepancies across the whole stack.
    /// Relative water scales use changes from the last accepted state.
    ///
    /// The controller charges both implicit stages and retains only accepted
    /// fine endpoint fluxes in each layer's water ledger. Interface transfers
    /// cancel in the whole-stack exterior account. Activity-range refusals
    /// shorten the proposal without clipping or first-order substitution.
    /// Faces, geometry and temperature stay fixed; local estimates do not
    /// certify global accuracy. Any late failure discards the complete duration.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_sdirk2_adaptive_to(
        &self,
        target_s: f64,
        faces: &[MoistureLayeredFace; 2],
        step: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            LayeredMoistureAdvanceReport<MoistureLayeredSdirk2Report>,
        ),
        MoistureAdvanceError<MoistureError>,
    > {
        let (state, report) = evolve_fields_with_method(
            self,
            target_s,
            time,
            TimeMethod::Sdirk2,
            |state, dt| {
                state
                    .advance_sdirk2(dt, faces, step, nonlinear, cx)
                    .map(|(state, report)| (state, LayeredStep(report)))
            },
            |state| Ok(observe(state)),
            |old, coarse, fine| {
                water_distance(
                    old.layers()
                        .iter()
                        .zip(coarse.layers())
                        .zip(fine.layers())
                        .map(|((a, b), c)| (a.slab().field(), b.slab().field(), c.slab().field())),
                    time,
                    cx,
                )
            },
            cx,
        )?;
        Ok((state, finish(report, time.water_balance_tolerance_kg, cx)?))
    }
}

//! Accuracy control and actual per-layer water accounts for curved composites.

use super::*;
use crate::moisture::radial::layered::MoistureLayeredRadial;

pub(super) fn observe_radial(body: &MoistureLayeredRadial) -> Snapshot<'_> {
    let fields = body
        .layers()
        .iter()
        .enumerate()
        .map(|(i, layer)| (i as u64, layer.radial().field()))
        .collect();
    let mut offset = 0;
    let activities = body
        .layers()
        .iter()
        .enumerate()
        .map(|(i, layer)| {
            let end = offset + layer.radial().moisture().len();
            let values = &body.activity()[offset..end];
            offset = end;
            (i as u64, values)
        })
        .collect();
    Snapshot {
        fields: Specimens::Many(fields),
        field_activities: Some(Specimens::Many(activities)),
        reservoir_water_kg: None,
        activity: None,
    }
}

pub(super) fn radial_water_distance<E>(
    old: &MoistureLayeredRadial,
    coarse: &MoistureLayeredRadial,
    fine: &MoistureLayeredRadial,
    time: MoistureTimeControls,
    cx: &Cx<'_>,
) -> Result<f64, MoistureAdvanceError<E>> {
    water_distance(
        old.layers()
            .iter()
            .zip(coarse.layers())
            .zip(fine.layers())
            .map(|((a, b), c)| (a.radial().field(), b.radial().field(), c.radial().field())),
        time,
        cx,
    )
}

impl MoistureLayeredRadial {
    /// Complete a coated cylinder/sphere duration with original nonlinear BE
    /// endpoints. Each layer retains its sources, unequal shell masses and
    /// physical interface laws. Error comparisons resolve every local content
    /// and activity plus the L1 sum of physical shell-water discrepancies.
    ///
    /// `reservoir_water_absolute_tolerance_kg` supplies that whole-body water
    /// error scale even with sealed exteriors. Only accepted fine interface
    /// transfers enter each layer's fixed original-to-final balance. The whole
    /// body reports exterior exchange separately. Any later failure or
    /// cancellation publishes no accepted prefix.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive_to(
        &self,
        target_s: f64,
        faces: &[MoistureLayeredFace; 2],
        step: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, LayeredMoistureAdvanceReport), MoistureAdvanceError<MoistureError>> {
        let (state, report) = evolve_fields_with_method(
            self,
            target_s,
            time,
            TimeMethod::FirstOrder,
            |state, dt| {
                state
                    .advance(dt, faces, step, nonlinear, cx)
                    .map(|(state, report)| (state, LayeredStep(report)))
            },
            |state| Ok(observe_radial(state)),
            |old, coarse, fine| radial_water_distance(old, coarse, fine, time, cx),
            cx,
        )?;
        Ok((state, finish(report, time.water_balance_tolerance_kg, cx)?))
    }

    /// Complete a heterogeneous radial duration with native SDIRK2 and the
    /// original second-order step-doubling controller. Actual first-stage face
    /// loads remain inside each native endpoint; they are not counted again in
    /// the accepted exterior or per-layer ledgers.
    ///
    /// Both stages spend the work cap even on an early failure. The comparison
    /// scales from [`Self::advance_adaptive_to`] apply with cubic local-error
    /// scaling. An inadmissible activity stage is retried at shorter duration,
    /// without clipping or a different transport model. Geometry, temperature
    /// and face conditions remain held; cut calls at physical events.
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
            |state| Ok(observe_radial(state)),
            |old, coarse, fine| radial_water_distance(old, coarse, fine, time, cx),
            cx,
        )?;
        Ok((state, finish(report, time.water_balance_tolerance_kg, cx)?))
    }
}

impl FiniteMoistureAtmosphere {
    /// Evolve one finite atmosphere and all original layers of a coated radial
    /// body to the same target. Every BE trial solves actual gas humidity and
    /// native heterogeneous diffusion together. Temporal control compares gas
    /// water/RH, each source activity/content and local dry-mass water error.
    ///
    /// The original controller checks both shell/interface transport and total
    /// vapor-plus-solid water. Final checks also enclose total and per-carrier
    /// gas mass over the entire duration. Any late failure discards every body
    /// and gas candidate. This adds neither a thermal solve nor condensation.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_layered_radial_adaptive_to(
        &self,
        moisture: &MoistureLayeredRadial,
        target_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        physical: AtmosphereMoistureSolve,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            MoistureLayeredRadial,
            LayeredMoistureAdvanceReport<AtmosphereMoistureStep<MoistureLayeredReport>>,
        ),
        MoistureAdvanceError<AtmosphereMoistureError>,
    > {
        checkpoint(cx)?;
        let initial = (self.clone(), moisture.clone());
        let ((atmosphere, field), report) = evolve_fields_with_method(
            &initial,
            target_s,
            time,
            TimeMethod::FirstOrder,
            |state, dt| {
                let (air, field, step) = state
                    .0
                    .advance_layered_radial(
                        &state.1, dt, exposed, transport, nonlinear, physical, cx,
                    )?
                    .into_parts();
                Ok(((air, field), LayeredStep(step)))
            },
            |state| {
                let mut snapshot = observe_radial(&state.1);
                snapshot.reservoir_water_kg = Some(state.0.vapor_mass_kg()?);
                snapshot.activity = Some(state.0.relative_humidity()?);
                Ok(snapshot)
            },
            |old, coarse, fine| radial_water_distance(&old.1, &coarse.1, &fine.1, time, cx),
            cx,
        )?;
        let mut report = finish(report, time.water_balance_tolerance_kg, cx)?;
        finite::check_vapor(
            self,
            &atmosphere,
            report.summary.water_change_kg,
            time.water_balance_tolerance_kg,
            &mut report.summary,
            cx,
        )?;
        Ok((atmosphere, field, report))
    }
}

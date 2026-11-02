//! Original heterogeneous layers and finite vapor under the shared time driver.
//!
//! Each original layer remains an independent inventory view with its own dry
//! cell mass. Its actual interface transfer enters adjacent layer balances with
//! opposite signs. Only the two exterior transfers enter the stack account.

use super::*;
use crate::moisture::{
    atmosphere::{
        AtmosphereMoistureError, AtmosphereMoistureSolve, AtmosphereMoistureStep,
        FiniteMoistureAtmosphere,
    },
    layered::{MoistureLayeredFace, MoistureLayeredReport, MoistureLayeredSlab},
};

/// Whole-stack accuracy and water accounts, plus each original material layer.
#[derive(Clone, Debug)]
pub struct LayeredMoistureAdvanceReport<Step = MoistureLayeredReport> {
    /// Face slots are the actual left/right exteriors, excluding internal joins.
    pub summary: MoistureAdvanceReport<Step, 2>,
    /// Physical left-to-right layer order; `id` is the zero-based layer index.
    /// Each layer's boundary account includes its actual material joins. The
    /// error ratio includes both dry-basis moisture and local source activity.
    pub layers: Vec<MoistureSpecimenAdvanceReport<2>>,
}

trait LayeredTransport {
    fn transport(&self) -> &MoistureLayeredReport;
}
impl LayeredTransport for MoistureLayeredReport {
    fn transport(&self) -> &MoistureLayeredReport {
        self
    }
}
impl LayeredTransport for AtmosphereMoistureStep<MoistureLayeredReport> {
    fn transport(&self) -> &MoistureLayeredReport {
        &self.transport
    }
}
struct LayeredStep<R>(R);
impl<R: LayeredTransport> ReportTransfers<2> for LayeredStep<R> {
    fn transfers(&self) -> Specimens<(u64, [f64; 2])> {
        let step = self.0.transport();
        let last = step.interfaces.len();
        Specimens::Many(
            (0..=last)
                .map(|layer| {
                    let left = if layer == 0 {
                        step.boundary_water_kg[0]
                    } else {
                        step.interfaces[layer - 1].water_left_to_right_kg
                    };
                    let right = if layer == last {
                        step.boundary_water_kg[1]
                    } else {
                        -step.interfaces[layer].water_left_to_right_kg
                    };
                    (layer as u64, [left, right])
                })
                .collect(),
        )
    }
}

fn observe(stack: &MoistureLayeredSlab) -> Snapshot<'_> {
    let fields = stack
        .layers()
        .iter()
        .enumerate()
        .map(|(i, layer)| (i as u64, layer.slab().field()))
        .collect();
    let mut offset = 0;
    let activities = stack
        .layers()
        .iter()
        .enumerate()
        .map(|(i, layer)| {
            let end = offset + layer.slab().moisture().len();
            let values = &stack.activity()[offset..end];
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

fn finish<R, E>(
    report: SharedMoistureAdvanceReport<LayeredStep<R>, 2>,
    tolerance: f64,
    cx: &Cx<'_>,
) -> Result<LayeredMoistureAdvanceReport<R>, MoistureAdvanceError<E>> {
    let mut summary = report.summary.map_step(|step| step.0);
    let layers = report.specimens;
    // The generic driver's per-layer faces include internal exchanges. Keep
    // those layer accounts, but expose only physical exterior stack transfers.
    summary.boundary_water_kg = [
        layers
            .first()
            .expect("admitted nonempty stack")
            .boundary_water_kg[0],
        layers
            .last()
            .expect("admitted nonempty stack")
            .boundary_water_kg[1],
    ];
    let mut transferred = WaterSum::default();
    for value in summary.boundary_water_kg {
        transferred.add(value);
    }
    let residual = summary.water_change_kg - transferred.total();
    if !residual.is_finite() || residual.abs() > tolerance {
        return Err(MoistureAdvanceError::WaterBalance {
            residual_kg: residual,
            tolerance_kg: tolerance,
        });
    }
    summary.transport_balance_residual_kg = residual;
    if summary.reservoir_water_kg.is_none() {
        summary.water_balance_residual_kg = residual;
    }
    checkpoint(cx)?;
    Ok(LayeredMoistureAdvanceReport { summary, layers })
}

impl MoistureLayeredSlab {
    /// Advance the complete heterogeneous stack to an exact target time.
    /// Original sorption, Kirchhoff diffusivity and contact laws are evaluated
    /// by the native layered owner on every full/half probe. Per-cell moisture
    /// and activity errors have independent physical scales. Accepted interface
    /// transfers audit each layer as well as the exterior whole-stack balance.
    pub fn advance_adaptive_to(
        &self,
        target_s: f64,
        faces: &[MoistureLayeredFace; 2],
        step: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        time: MoistureTimeControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, LayeredMoistureAdvanceReport), MoistureAdvanceError<MoistureError>> {
        let (state, report) = evolve_fields(
            self,
            target_s,
            time,
            |state, dt| {
                state
                    .advance(dt, faces, step, nonlinear, cx)
                    .map(|(state, step)| (state, LayeredStep(step)))
            },
            |state| Ok(observe(state)),
            cx,
        )?;
        Ok((state, finish(report, time.water_balance_tolerance_kg, cx)?))
    }
}

impl FiniteMoistureAtmosphere {
    /// Advance one original finite gas inventory and every heterogeneous layer
    /// together. Source/temperature/contact laws stay fixed. RH, vapor mass,
    /// local activity and dry-basis moisture all enter temporal accuracy; actual
    /// carrier and total mass are independently audited over the whole window.
    /// Any failure, including a late source refusal or cancellation, returns no
    /// partially advanced gas, layers, clock or source receipts.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_layered_adaptive_to(
        &self,
        moisture: &MoistureLayeredSlab,
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
            MoistureLayeredSlab,
            LayeredMoistureAdvanceReport<AtmosphereMoistureStep<MoistureLayeredReport>>,
        ),
        MoistureAdvanceError<AtmosphereMoistureError>,
    > {
        checkpoint(cx)?;
        let initial = (self.clone(), moisture.clone());
        let ((atmosphere, field), report) = evolve_fields(
            &initial,
            target_s,
            time,
            |state, dt| {
                let (air, field, step) = state
                    .0
                    .advance_layered(&state.1, dt, exposed, transport, nonlinear, physical, cx)?
                    .into_parts();
                Ok(((air, field), LayeredStep(step)))
            },
            |state| {
                let mut snapshot = observe(&state.1);
                snapshot.reservoir_water_kg = Some(state.0.vapor_mass_kg()?);
                snapshot.activity = Some(state.0.relative_humidity()?);
                Ok(snapshot)
            },
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

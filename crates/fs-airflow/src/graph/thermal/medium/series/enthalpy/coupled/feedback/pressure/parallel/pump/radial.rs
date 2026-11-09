//! An original sourced pump curve in the native transient radial graph.
//!
//! The pump's original material query and equipment object stay shared. Head
//! is converted with that inlet density and original gravity, while the exact
//! piecewise-linear tangent enters the native pressure Schur complement. No
//! affinity law, reversed operating quadrant, efficiency or thermal loss is
//! inferred. Optional original drive-curve refusals remain active.

use std::sync::Arc;

use fs_material::liquid::ResolvedLiquidState;
use fs_qty::{Pressure, VolumetricFlowRate};

use super::{PumpHeadCurve, PumpOperatingError, PumpSupplyState};
use crate::graph::rheology::transient::head::power::FlowHeadPowerReceipt;
use crate::graph::rheology::transient::head::{FlowHeadLaw, FlowHeadResponse, FlowHeadSource};

/// Source-preserving adapter for a fixed prepared equipment/material state.
/// The radial graph requires its constant reference density to equal this
/// original inlet density. Regional fixed-mass thermal transport retains that
/// same restriction: refreshing this material is not an expansion/EOS model.
#[derive(Clone, Debug)]
pub struct RadialPumpHead {
    curve: Arc<PumpHeadCurve>,
    inlet: Arc<ResolvedLiquidState>,
    support: [f64; 2],
}

impl RadialPumpHead {
    /// Keep the exact original owners, admitting their inlet-viscosity support.
    /// When drive data are bound, the flow span is the intersection of the two
    /// authored source spans. Actual drive-power/material feasibility is still
    /// checked by the original pump at every queried flow, including brackets.
    pub fn try_new(
        curve: Arc<PumpHeadCurve>,
        inlet: Arc<ResolvedLiquidState>,
    ) -> Result<Self, PumpOperatingError> {
        curve.admit_inlet(&inlet)?;
        let points = curve.points();
        let mut support = [
            points[0].flow.value(),
            points[points.len() - 1].flow.value(),
        ];
        if let Some(drive) = curve.drive_curve() {
            let points = drive.points();
            support[0] = support[0].max(points[0].flow.value());
            support[1] = support[1].min(points[points.len() - 1].flow.value());
        }
        if support[0] >= support[1] {
            return Err(PumpOperatingError::Invalid(
                "original pump and drive require an overlapping flow interval",
            ));
        }
        Ok(Self {
            curve,
            inlet,
            support,
        })
    }

    /// The same equipment object, including its provenance, speed and drive data.
    pub fn curve(&self) -> &PumpHeadCurve {
        &self.curve
    }
    /// The same material query and retained original property-use receipts.
    pub fn inlet_material(&self) -> &ResolvedLiquidState {
        &self.inlet
    }
    /// Share this prepared source directly with native or regional drive data.
    pub fn source(self: &Arc<Self>) -> FlowHeadSource {
        FlowHeadSource::new(self.clone())
    }

    /// Original hydraulic and optional shaft/electrical supply at the actual
    /// solved flow. Native accepted head responses already retain this same
    /// evaluation's power state; consumers need not query again to route losses.
    pub fn supply(&self, flow: VolumetricFlowRate) -> Result<PumpSupplyState, PumpOperatingError> {
        self.curve.evaluate(&self.inlet, flow)
    }
}

impl FlowHeadLaw for RadialPumpHead {
    fn flow_interval_m3_s(&self) -> [f64; 2] {
        self.support
    }

    fn reference_density_kg_m3(&self) -> Option<f64> {
        Some(self.inlet.density_kg_m3())
    }

    fn evaluate(
        &self,
        flow: VolumetricFlowRate,
    ) -> Result<FlowHeadResponse, Box<dyn std::error::Error + Send + Sync>> {
        let supply = self.supply(flow)?;
        let points = self.curve.points();
        // At an interior knot use the right derivative; at the last endpoint
        // use the left derivative. No averaged or smoothed curve is introduced.
        let upper = points
            .partition_point(|point| point.flow.value() <= flow.value())
            .min(points.len() - 1);
        let left = points[upper - 1];
        let right = points[upper];
        let slope = supply.inlet_density_kg_m3
            * self.curve.gravity_m_per_s2()
            * ((right.head.value() - left.head.value()) / (right.flow.value() - left.flow.value()));
        if !slope.is_finite() || slope > 0.0 {
            return Err(PumpOperatingError::Invalid(
                "unrepresentable original pump pressure tangent",
            )
            .into());
        }
        Ok(FlowHeadResponse {
            pressure: Pressure::new(supply.available_pressure_pa),
            pressure_slope_pa_s_m3: slope,
            power: supply.drive_power.map(|state| FlowHeadPowerReceipt {
                state,
                routing: self
                    .curve
                    .drive_curve()
                    .expect("original supply has drive data")
                    .routing(),
            }),
        })
    }
}

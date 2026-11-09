//! Fixed-speed equipment input curves, including nonzero shutoff consumption.
//!
//! Head/flow data determine hydraulic output, not input power. Shaft and
//! electrical inputs therefore come from separate authored measurements/models
//! at the same speed, over explicit density/viscosity domains. No constant
//! efficiency, affinity-law scaling, or motor loss is invented. Their
//! piecewise-linear interpolation must satisfy electrical >= shaft >= hydraulic
//! at every actual trial, including between the original samples.

use super::{PumpOperatingError, SourceProvenance};
use fs_qty::{Power, VolumetricFlowRate};

// The source evaluator and native head owner use the same original value
// types. Preserve the established public pump API without a second schema.
pub use crate::graph::rheology::transient::head::power::{PumpLossRouting, PumpPowerState};

/// One equipment point at the head curve's fixed speed and pump-inlet flow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PumpPowerPoint {
    /// Nonnegative volume flow at the pump inlet.
    pub flow: VolumetricFlowRate,
    /// Mechanical power delivered by the drive to the pump shaft.
    pub shaft_input: Power,
    /// Electrical power consumed by the drive, including its losses.
    pub electrical_input: Power,
}

/// Independently sourced equipment data, with no extrapolation or scaling.
#[derive(Clone, Debug, PartialEq)]
pub struct PumpDriveCurve {
    points: Vec<PumpPowerPoint>,
    density_domain_kg_m3: (f64, f64),
    viscosity_domain_m2_s: (f64, f64),
    routing: PumpLossRouting,
    source: SourceProvenance,
}

impl PumpDriveCurve {
    /// Admit 2..=1024 increasing flows, nonnegative input powers, electrical
    /// input >= shaft input, positive inclusive fluid domains and explicit heat
    /// fractions in [0,1]. Input power may be nonzero at zero hydraulic flow.
    /// The source must justify use at the head curve's speed and across BOTH
    /// fluid domains. This is a nominal equipment model, not a measured guarantee.
    pub fn try_new(
        points: Vec<PumpPowerPoint>,
        density_domain_kg_m3: (f64, f64),
        viscosity_domain_m2_s: (f64, f64),
        routing: PumpLossRouting,
        source: SourceProvenance,
    ) -> Result<Self, PumpOperatingError> {
        if !(2..=1024).contains(&points.len())
            || points.iter().any(|p| {
                [
                    p.flow.value(),
                    p.shaft_input.value(),
                    p.electrical_input.value(),
                ]
                .iter()
                .any(|v| !v.is_finite() || *v < 0.0)
                    || p.electrical_input.value() < p.shaft_input.value()
            })
            || points
                .windows(2)
                .any(|p| p[1].flow.value() <= p[0].flow.value())
            || !valid_domain(density_domain_kg_m3)
            || !valid_domain(viscosity_domain_m2_s)
            || ![routing.pump_to_reservoir, routing.drive_to_reservoir]
                .iter()
                .all(|v| v.is_finite() && (0.0..=1.0).contains(v))
            || source.citation.trim().is_empty()
            || source.identifier.trim().is_empty()
        {
            return Err(PumpOperatingError::Invalid(
                "invalid drive power curve, fluid domain, heat routing or source",
            ));
        }
        Ok(Self {
            points,
            density_domain_kg_m3,
            viscosity_domain_m2_s,
            routing,
            source,
        })
    }

    /// Original equipment samples in ascending flow order.
    pub fn points(&self) -> &[PumpPowerPoint] {
        &self.points
    }
    /// Provenance of the power measurements or model and its applicability.
    pub const fn source(&self) -> &SourceProvenance {
        &self.source
    }
    /// Declared thermal destination fractions of the two distinct losses.
    pub const fn routing(&self) -> PumpLossRouting {
        self.routing
    }
    /// Inclusive pump-inlet density support in kg/m³.
    pub const fn density_domain_kg_m3(&self) -> (f64, f64) {
        self.density_domain_kg_m3
    }
    /// Inclusive pump-inlet kinematic viscosity support in m²/s.
    pub const fn viscosity_domain_m2_s(&self) -> (f64, f64) {
        self.viscosity_domain_m2_s
    }

    // Reused by calibrated-speed preparation; hydraulic feasibility still
    // belongs to evaluate at the actual material and final commanded speed.
    pub(super) fn inputs_at_flow(&self, q: f64) -> Result<(f64, f64), PumpOperatingError> {
        in_domain(
            "drive inlet volume flow",
            q,
            (
                self.points[0].flow.value(),
                self.points[self.points.len() - 1].flow.value(),
            ),
        )?;
        let upper = self.points.partition_point(|p| p.flow.value() < q);
        let right = self.points[upper];
        Ok(if right.flow.value() == q {
            (right.shaft_input.value(), right.electrical_input.value())
        } else {
            let left = self.points[upper - 1];
            let fraction = (q - left.flow.value()) / (right.flow.value() - left.flow.value());
            (
                (1.0 - fraction) * left.shaft_input.value() + fraction * right.shaft_input.value(),
                (1.0 - fraction) * left.electrical_input.value()
                    + fraction * right.electrical_input.value(),
            )
        })
    }

    pub(super) fn evaluate(
        &self,
        flow: VolumetricFlowRate,
        density: f64,
        viscosity: f64,
        hydraulic_power_w: f64,
    ) -> Result<PumpPowerState, PumpOperatingError> {
        in_domain("drive inlet density", density, self.density_domain_kg_m3)?;
        in_domain(
            "drive inlet kinematic viscosity",
            viscosity,
            self.viscosity_domain_m2_s,
        )?;
        let (shaft, electrical) = self.inputs_at_flow(flow.value())?;
        // A head/power pair can disagree only BETWEEN samples: hydraulic power
        // is quadratic in flow on a linear head segment. Never clip that deficit
        // or silently reinterpret a physically impossible equipment curve.
        if !hydraulic_power_w.is_finite()
            || hydraulic_power_w < 0.0
            || shaft < hydraulic_power_w
            || electrical < shaft
        {
            return Err(PumpOperatingError::Invalid(
                "drive input cannot supply the actual hydraulic output",
            ));
        }
        let pump_loss = shaft - hydraulic_power_w;
        let drive_loss = electrical - shaft;
        let reservoir = self.routing.pump_to_reservoir * pump_loss
            + self.routing.drive_to_reservoir * drive_loss;
        let surroundings = (1.0 - self.routing.pump_to_reservoir) * pump_loss
            + (1.0 - self.routing.drive_to_reservoir) * drive_loss;
        if [
            shaft,
            electrical,
            pump_loss,
            drive_loss,
            reservoir,
            surroundings,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(PumpOperatingError::Invalid(
                "drive power or heat partition is unrepresentable",
            ));
        }
        Ok(PumpPowerState {
            shaft_input_w: shaft,
            electrical_input_w: electrical,
            pump_loss_w: pump_loss,
            drive_loss_w: drive_loss,
            reservoir_heat_w: reservoir,
            surroundings_heat_w: surroundings,
        })
    }
}

fn valid_domain((lo, hi): (f64, f64)) -> bool {
    lo.is_finite() && hi.is_finite() && lo > 0.0 && hi >= lo
}
fn in_domain(
    field: &'static str,
    value: f64,
    (lower, upper): (f64, f64),
) -> Result<(), PumpOperatingError> {
    if value.is_finite() && value >= lower && value <= upper {
        Ok(())
    } else {
        Err(PumpOperatingError::OutsideCurve {
            field,
            value,
            lower,
            upper,
        })
    }
}

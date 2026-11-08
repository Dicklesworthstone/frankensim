//! Source-preserving turbulent pipe branches for the native passive graph.
//!
//! The nodal solver brackets through zero flow. Outside the authored turbulent
//! Reynolds band it uses a private monotone numerical continuation, solely to
//! locate roots. Every public graph/fan result re-evaluates the original pipe;
//! continuation is never an admissible laminar or transition model.

use fs_material::fluid_transport::FluidTransportState;
use fs_qty::{MassFlowRate, VolumetricFlowRate};

use crate::graph::thermal::medium::series::hydraulic::{
    MaterialFlowResistance, ResistanceError, ResistanceEvaluation,
};
use crate::{LossElement, LossResistance, SourceProvenance, ToleranceBasis};

/// One original Colebrook pipe and its coherent frozen material snapshot.
/// Momentum density is explicit for the existing bounded Boussinesq owner;
/// no material property, roughness, or friction coefficient is fitted here.
#[derive(Debug, Clone, PartialEq)]
pub struct TurbulentPipeLoss {
    name: String,
    law: MaterialFlowResistance,
    transport: FluidTransportState,
    momentum_density_kg_m3: f64,
    uncertainty_rel: f64,
    source: SourceProvenance,
    tolerance_basis: ToleranceBasis,
    lower_flow_m3_s: f64,
    upper_flow_m3_s: f64,
    lower: ResistanceEvaluation,
    upper: ResistanceEvaluation,
}

impl TurbulentPipeLoss {
    /// Bind a developed turbulent pipe to one original material snapshot.
    /// The source describes its geometric law and developed-flow assumption.
    /// Source density bounds remain the caller's explicit Boussinesq policy.
    ///
    /// # Errors
    /// A non-Colebrook law, missing source, invalid density/uncertainty, or
    /// unrepresentable Reynolds endpoints, pressure loss, power or tangent.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        name: impl Into<String>,
        law: MaterialFlowResistance,
        transport: FluidTransportState,
        momentum_density_kg_m3: f64,
        uncertainty_rel: f64,
        source: SourceProvenance,
        tolerance_basis: ToleranceBasis,
    ) -> Result<Self, ResistanceError> {
        let name = name.into();
        if name.trim().is_empty() || source.citation.trim().is_empty()
            || source.identifier.trim().is_empty()
            || !uncertainty_rel.is_finite() || !(0.0..1.0).contains(&uncertainty_rel)
            || !momentum_density_kg_m3.is_finite() || momentum_density_kg_m3 <= 0.0
            || law.relative_sand_roughness().is_none()
        {
            return Err(ResistanceError::Invalid("source-backed turbulent pipe and positive momentum density required"));
        }
        // Round the two numerical brackets inward, never widen source support.
        // The original evaluation below remains the authority on every endpoint.
        let (minimum, maximum) = law.reynolds_range();
        let endpoint = |re: f64, lower: bool| {
            let mut q = ((re * transport.dynamic_viscosity_pa_s())
                * (law.flow_area().value() / law.reynolds_length().value())) / momentum_density_kg_m3;
            for _ in 0..16 {
                let mass = q * momentum_density_kg_m3;
                match law.evaluate_snapshot(transport, MassFlowRate::new(mass), momentum_density_kg_m3) {
                    Ok(value) => return Ok((q, value)),
                    Err(ResistanceError::Reynolds { value, .. })
                        if (lower && value < minimum) || (!lower && value > maximum) => {
                            q = if lower { q.next_up() } else { q.next_down() };
                        }
                    Err(error) => return Err(error),
                }
            }
            Err(ResistanceError::Invalid("turbulent Reynolds endpoint is unrepresentable"))
        };
        let (lower_flow_m3_s, lower) = endpoint(minimum, true)?;
        let (upper_flow_m3_s, upper) = endpoint(maximum, false)?;
        if lower_flow_m3_s >= upper_flow_m3_s
            || lower.pressure_drop_pa >= upper.pressure_drop_pa
            || !(upper.dpressure_dmass_flow * momentum_density_kg_m3).is_finite()
            || upper.dpressure_dmass_flow * momentum_density_kg_m3 <= 0.0
        {
            return Err(ResistanceError::Invalid("turbulent numerical bracket is unrepresentable"));
        }
        Ok(Self { name, law, transport, momentum_density_kg_m3, uncertainty_rel,
            source, tolerance_basis, lower_flow_m3_s, upper_flow_m3_s, lower, upper })
    }

    /// Original geometry, equivalent sand roughness and authored Reynolds band.
    #[must_use]
    pub const fn law(&self) -> MaterialFlowResistance { self.law }

    /// Coherent material properties used by this fixed-state branch.
    #[must_use]
    pub const fn transport(&self) -> FluidTransportState { self.transport }

    /// Conserved mass/volume conversion and momentum density.
    #[must_use]
    pub const fn momentum_density_kg_m3(&self) -> f64 { self.momentum_density_kg_m3 }

    /// Original geometric-law provenance; uncertainty is retained, not propagated.
    #[must_use]
    pub const fn source(&self) -> &SourceProvenance { &self.source }

    /// Evaluate the ORIGINAL law, including its Reynolds support and power.
    /// Zero/transition flow refuses; the graph's numerical continuation is private.
    pub fn evaluate(&self, flow: VolumetricFlowRate) -> Result<ResistanceEvaluation, ResistanceError> {
        let q = flow.value();
        let mass = q * self.momentum_density_kg_m3;
        if !q.is_finite() || !mass.is_finite() || (q != 0.0 && mass == 0.0) {
            return Err(ResistanceError::Invalid("turbulent graph mass flow is unrepresentable"));
        }
        self.law.evaluate_snapshot(self.transport, MassFlowRate::new(mass), self.momentum_density_kg_m3)
    }

    pub(super) fn numerical_flow(&self, pressure: f64) -> Result<f64, ResistanceError> {
        let dp = pressure.abs();
        if !dp.is_finite() {
            return Err(ResistanceError::Invalid("nonfinite turbulent branch pressure"));
        }
        let q = if dp == 0.0 {
            0.0
        } else if dp < self.lower.pressure_drop_pa {
            // Positive secant to the origin, used only inside root brackets.
            self.lower_flow_m3_s * (dp / self.lower.pressure_drop_pa)
        } else if dp > self.upper.pressure_drop_pa {
            // Positive endpoint tangent avoids a false high-flow saturation.
            self.upper_flow_m3_s + ((dp - self.upper.pressure_drop_pa)
                / self.upper.dpressure_dmass_flow) / self.momentum_density_kg_m3
        } else {
            // Suppress only inverse-rounding excursions at the two admitted
            // endpoints. Original support and pressure closure are still checked.
            self.law.turbulent_trial_volume_flow(self.transport, self.momentum_density_kg_m3, dp)?
                .clamp(self.lower_flow_m3_s, self.upper_flow_m3_s)
        };
        if !q.is_finite() || (dp != 0.0 && q <= 0.0) {
            return Err(ResistanceError::Invalid("turbulent graph trial flow is unrepresentable"));
        }
        Ok(q.copysign(pressure))
    }

    pub(super) fn zero_quadratic_metadata(&self) -> LossElement {
        LossElement {
            name: self.name.clone(), resistance: LossResistance::new(0.0),
            uncertainty_rel: self.uncertainty_rel, source: self.source.clone(),
            tolerance_basis: self.tolerance_basis, regime_validity: None,
        }
    }
}

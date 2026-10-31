//! Darcy and calibrated Forchheimer drag in the existing passive pressure graph.
//!
//! This graph conserves volume at one frozen, low-Mach gas state. A porous
//! element uses the material owner's intrinsic permeability and actual gas
//! viscosity; it is not a compressible mass-flow or membrane-diffusion model.
//! Source velocity bounds and the caller's small-pressure-drop allowance are
//! checked on the accepted solution, including a fan-driven solution.
//! Graph pressures are gauge pressures relative to the common frozen gas
//! state's absolute pressure. Each face must remain within both the material's
//! absolute-pressure support and the allowed density excursion from that state.

use fs_material::{gas::GasState, porous::permeability::DarcyGasBarrier};
use fs_qty::Qty;

use super::{GraphBranch, GraphError};
use crate::{LossElement, LossResistance, SourceProvenance, ToleranceBasis};

mod forchheimer;

/// Linear hydraulic resistance [Pa s/m³], distinct from quadratic resistance.
pub type LinearLossResistance = Qty<-4, 1, -1, 0, 0, 0>;

/// The retained source material and operating state behind a Darcy coefficient.
#[derive(Debug, Clone, PartialEq)]
pub struct DarcyLossBinding {
    /// Original layers, orientation, and source support.
    pub barrier: DarcyGasBarrier,
    /// Frozen common gas state; its pressure is the absolute pressure at graph
    /// gauge zero. A face's absolute pressure is this pressure plus its node pressure.
    pub gas: GasState,
    /// Exposed superficial-flow area [m²].
    pub area_m2: f64,
    /// Explicit allowance for both `|dp| / p_reference` and each face's
    /// `|p_gauge| / p_reference` in this constant-density model.
    pub maximum_relative_pressure_drop: f64,
}

/// A passive linear pressure loss with retained physical authority.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearLossElement {
    /// Stable component identity.
    pub name: String,
    /// Nominal coefficient in `dp = R_linear Q` [Pa s/m³].
    pub resistance: LinearLossResistance,
    /// Declared coefficient uncertainty, retained but not propagated by the graph.
    pub uncertainty_rel: f64,
    /// Authority for the coefficient and its uncertainty.
    pub source: SourceProvenance,
    /// Authority classification for the uncertainty.
    pub tolerance_basis: ToleranceBasis,
    material: Option<DarcyLossBinding>,
    // Some requires this exact quadratic coefficient on the same graph edge.
    // None preserves independently authored linear/quadratic compositions.
    quadratic_binding: Option<LossResistance>,
}

impl LinearLossElement {
    /// Construct a source-declared linear hydraulic resistance.
    pub fn new(
        name: impl Into<String>,
        resistance: LinearLossResistance,
        uncertainty_rel: f64,
        source: SourceProvenance,
        tolerance_basis: ToleranceBasis,
    ) -> Result<Self, GraphError> {
        let element = Self {
            name: name.into(), resistance, uncertainty_rel, source, tolerance_basis,
            material: None, quadratic_binding: None,
        };
        element.validate().map_err(GraphError::InvalidInput)?;
        Ok(element)
    }

    /// Resolve the shared material law into frozen-property Darcy resistance.
    ///
    /// The material owner checks temperature/pressure support, gas coherence,
    /// geometry, and absence of a slip correction. A sealed barrier is not an
    /// epsilon-conductance edge: omit that path from the pressure graph instead.
    /// No inertial coefficient is inferred from permeability.
    #[allow(clippy::too_many_arguments)]
    pub fn from_darcy(
        name: impl Into<String>,
        barrier: DarcyGasBarrier,
        gas: GasState,
        area_m2: f64,
        maximum_relative_pressure_drop: f64,
        uncertainty_rel: f64,
        source: SourceProvenance,
        tolerance_basis: ToleranceBasis,
    ) -> Result<Self, GraphError> {
        if !(maximum_relative_pressure_drop.is_finite()
            && maximum_relative_pressure_drop > 0.0
            && maximum_relative_pressure_drop < 1.0)
        {
            return Err(GraphError::InvalidInput("Darcy relative pressure allowance must be in (0, 1)"));
        }
        let resistance = barrier.linear_volume_resistance(&gas, area_m2)
            .map_err(GraphError::Material)?
            .ok_or(GraphError::InvalidInput("a sealed Darcy path must be omitted from the graph"))?;
        let mut element = Self::new(name, LinearLossResistance::new(resistance),
            uncertainty_rel, source, tolerance_basis)?;
        element.material = Some(DarcyLossBinding {
            barrier, gas, area_m2, maximum_relative_pressure_drop,
        });
        Ok(element)
    }

    /// Complete original material binding, when constructed from permeability.
    #[must_use]
    pub fn material(&self) -> Option<&DarcyLossBinding> { self.material.as_ref() }

    pub(super) fn validate(&self) -> Result<(), &'static str> {
        if self.name.trim().is_empty() {
            return Err("linear loss name must be nonblank");
        }
        if !(self.resistance.value().is_finite() && self.resistance.value() > 0.0) {
            return Err("linear resistance must be finite and positive");
        }
        if !(self.uncertainty_rel.is_finite() && (0.0..1.0).contains(&self.uncertainty_rel)) {
            return Err("linear resistance uncertainty must be in [0, 1)");
        }
        if self.source.citation.trim().is_empty() || self.source.identifier.trim().is_empty() {
            return Err("linear loss needs an explicit source");
        }
        if let Some(binding) = &self.material {
            let expected = if let Some(quadratic) = self.quadratic_binding {
                let coefficients = binding.barrier
                    .volume_resistance_coefficients(&binding.gas, binding.area_m2)
                    .map_err(|_| "porous material binding is outside its source support")?
                    .ok_or("a sealed porous path must be omitted from the graph")?;
                if coefficients.1 != quadratic.value() {
                    return Err("porous quadratic resistance no longer matches its material");
                }
                coefficients.0
            } else {
                binding.barrier.linear_volume_resistance(&binding.gas, binding.area_m2)
                    .map_err(|_| "Darcy material binding is outside its source support")?
                    .ok_or("a sealed Darcy path must be omitted from the graph")?
            };
            if self.resistance.value() != expected {
                return Err("Darcy resistance no longer matches its retained material binding");
            }
        }
        Ok(())
    }

    pub(super) fn validate_quadratic(&self, actual: f64) -> Result<(), &'static str> {
        if self.quadratic_binding.is_some_and(|expected| expected.value() != actual) {
            return Err(
                "porous edge must retain its original calibrated quadratic resistance",
            );
        }
        Ok(())
    }

    pub(super) fn admit_flow(
        &self, from_gauge_pressure: f64, to_gauge_pressure: f64, flow: f64,
    ) -> Result<(), &'static str> {
        let Some(binding) = &self.material else { return Ok(()); };
        let pressure_drop = from_gauge_pressure - to_gauge_pressure;
        if !pressure_drop.is_finite() {
            return Err("Darcy pressure drop must be finite");
        }
        if pressure_drop.abs() / binding.gas.pressure > binding.maximum_relative_pressure_drop {
            return Err("Darcy pressure drop exceeds the frozen-density allowance");
        }
        if flow.abs() / binding.area_m2 > binding.barrier.support().maximum_superficial_speed_m_s {
            return Err("Darcy superficial velocity exceeds material source support");
        }
        let (minimum, maximum) = binding.barrier.support().pressure_pa;
        for gauge_pressure in [from_gauge_pressure, to_gauge_pressure] {
            let absolute_pressure = binding.gas.pressure + gauge_pressure;
            if !(absolute_pressure.is_finite() && absolute_pressure > 0.0) {
                return Err("Darcy endpoint absolute pressure must be finite and positive");
            }
            if absolute_pressure < minimum || absolute_pressure > maximum {
                return Err("Darcy endpoint pressure is outside material source support");
            }
            if gauge_pressure.abs() / binding.gas.pressure > binding.maximum_relative_pressure_drop {
                return Err("Darcy endpoint density excursion exceeds the frozen-state allowance");
            }
        }
        Ok(())
    }

    // Explicitly zero quadratic drag: internal graph metadata, never a fabricated
    // positive LossElement usable by the quadratic-only series/parallel owner.
    fn zero_quadratic_metadata(&self) -> LossElement {
        LossElement {
            name: self.name.clone(), resistance: LossResistance::new(0.0),
            uncertainty_rel: self.uncertainty_rel, source: self.source.clone(),
            tolerance_basis: self.tolerance_basis, regime_validity: None,
        }
    }
}

/// One graph edge with its actual passive pressure/volume-flow law.
#[derive(Debug, Clone, PartialEq)]
pub enum PassiveBranch {
    /// The existing `dp = R Q |Q|` law.
    Quadratic(GraphBranch),
    /// A pure linear law, including Darcy flow through a porous barrier.
    Linear {
        /// Positive-flow origin.
        from: usize,
        /// Positive-flow destination.
        to: usize,
        /// Source-backed linear drag.
        loss: LinearLossElement,
    },
    /// Independent linear and quadratic losses on the same flow path.
    ///
    /// A Forchheimer interpretation requires a source-backed quadratic
    /// coefficient; the permeability adapter does not invent one.
    LinearQuadratic {
        /// Existing quadratic component and endpoints.
        branch: GraphBranch,
        /// Additional linear component on that path.
        linear_loss: LinearLossElement,
    },
}

impl From<GraphBranch> for PassiveBranch {
    fn from(branch: GraphBranch) -> Self { Self::Quadratic(branch) }
}

impl PassiveBranch {
    pub(super) fn into_parts(self) -> (GraphBranch, Option<LinearLossElement>) {
        match self {
            Self::Quadratic(branch) => (branch, None),
            Self::Linear { from, to, loss } => (
                GraphBranch { from, to, loss: loss.zero_quadratic_metadata() }, Some(loss),
            ),
            Self::LinearQuadratic { branch, linear_loss } => (branch, Some(linear_loss)),
        }
    }
}

#[cfg(test)]
mod tests;
